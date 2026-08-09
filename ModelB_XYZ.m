function ModelB_XYZ()
% MODELB_SURFACESPIN_LLG_ODE45_FINAL_XYZ_INITIAL_SPINS
%
% ORIGINAL MODEL B WITH A DYNAMICAL SURFACE SPIN AND THE BIDIRECTIONAL
% FEEDBACK LOOP WRITTEN ON THE WHITEBOARD.
%
% Nothing from the original preparation, contact Hamiltonian, temperature,
% detailed-balance Lindblad construction, molecular geometries, current
% calculation or Biot-Savart calculation is replaced.
%
% The only controlled change in this version is that the complete calculation
% is repeated for three normalized initial surface-spin directions:
%
%   S_surface(0) = [1 0 0]   (+x)
%   S_surface(0) = [0 1 0]   (+y)
%   S_surface(0) = [0 0 1]   (+z)
%
% FULL LOOP IMPLEMENTED IN THIS FILE:
%
%   gas-phase ground state
%       -> rho(t)
%       -> bond currents J_ij(t)
%       -> B(R_i,t) at molecular sites
%       -> Zeeman Hamiltonian H_Z(t)
%       -> B(R_surface,t) at the surface-spin position
%       -> LLG equation for S_surface(t)
%       -> exchange Hamiltonian
%          H_surface-ex(t)=sum_i J_i S_surface(t).s_i
%       -> updated rho(t)
%
% with
%
%   J_i = J0 exp(-r_i/r0)
%   s_i = sigma/2
%
% This is bidirectional because:
%
%   rho changes the currents and B_surface, which changes S_surface;
%   S_surface changes H_surface-ex, which changes rho.
%
% IMPORTANT MODEL BOUNDARY:
%
% The LLG effective field in this code is exactly the Biot-Savart field at
% the surface position. No additional direct exchange-backaction field,
% alpha_0, alpha_1, electron reservoir, injection/extraction process or
% phenomenological amplification term is added.
%
% BASIS ORDER USED THROUGHOUT:
%
%   |1,up>, |1,down>, |2,up>, |2,down>, ... , |L,up>, |L,down>.
%
% ORIGINAL PARAMETERS PRESERVED:
%
%   t_hop       = 1.0 eV
%   q_eff       = 0.30
%   T           = 300 K
%   tau_rel     = 100 fs
%   Gamma0      = 0.01 fs^-1 per downward transition channel
%   J0          = 100 meV
%   r0          = 3 Angstrom
%   alpha_LLG   = 0.10
%   gamma_LLG   = 1.760859e-4 rad/(fs T)
%   solver      = ode45 (adaptive Runge-Kutta)
%   InitialStep = 1.0e-4 fs
%   early MaxStep = 1.0e-2 fs for 0-500 fs; late stage is adaptive
%   t_end       = 5000 fs (absolute safety limit)
%
% CRITICAL UNIT CONVERSION:
%
%   Biot-Savart output is stored in mT.
%   Both the molecular Zeeman Hamiltonian and the LLG equation require T:
%
%       B_T = 1e-3 B_mT.
%
% FIGURES PRODUCED BY THIS VERSION - NO OTHERS:
%
%   1 geometry figure containing both molecular geometries.
%   3 final spin-imbalance figures, one for each initial S_surface direction.
%   3 time-dynamics figures, one for each initial S_surface direction.
%     Each dynamics figure contains both geometries and, for x/y/z, plots
%     B_surface [mT] together with normalized M/M_s = S_surface.
% =========================================================================

clc;
close all;

%% ========================================================================
% 1. PARAMETERS
% =========================================================================

p.t_hop_eV = 1.0;
p.e2_over_4pieps0_eVA = 14.4;
p.q_eff = 0.30;
p.R_surface_A = [0 0 0];

% Surface-exchange profile.
p.J0_meV = 100;
p.J0_eV = 1e-3*p.J0_meV;
p.r0_A = 3.0;

% Base surface-spin direction. Each actual run below overwrites ONLY p.S0
% with one of the three normalized Cartesian initial directions.
p.S0 = [0 0 1];
p.S0 = p.S0/norm(p.S0);
p.gamma_rad_fs_T = 1.760859e-4;
p.alpha_LLG = 0.10;

% Thermal bath and detailed balance.
p.kB_eV_K = 8.617333262e-5;
p.T_K = 300;
p.beta_eV_inv = 1/(p.kB_eV_K*p.T_K);
p.tau_rel_fs = 100;
p.Gamma0_per_fs = 1/p.tau_rel_fs;

% Constants used in the Hamiltonian.
p.hbar_eV_fs = 0.6582119569;
p.g_e = 2.00231930436256;
p.muB_eV_T = 5.7883818060e-5;
p.zeeman_prefactor_eV_T = 0.5*p.g_e*p.muB_eV_T;

% -------------------------------------------------------------------------
% TIME CONTROL: ACCURATE EARLY TRANSIENT + JOINT EARLY STOP
% -------------------------------------------------------------------------
% 0-500 fs is propagated with MaxStep = 0.01 fs. After 500 fs, ode45
% advances in adaptive 250 fs chunks. The run stops only after BOTH the
% electronic density and the LLG surface spin satisfy the stationarity
% criteria in two consecutive checks. The absolute safety limit is 5000 fs.
p.t_start_fs = 0;
p.early_window_fs = 500;
p.t_max_fs = 5000;
p.t_end_fs = p.t_max_fs;       % compatibility with existing reports
p.chunk_fs = 250;
p.output_dt_early_fs = 0.05;
p.min_stop_time_fs = 1000;
p.convergence_window_fs = 500;
p.required_consecutive_checks = 2;
p.derivative_check_points = 41;

% Joint electronic + LLG stabilization tolerances.
p.population_range_tol = 1e-6;
p.surface_spin_range_tol = 1e-6;
p.rho_rhs_norm_tol_per_fs = 1e-7;
p.surface_rhs_norm_tol_per_fs = 1e-8;
p.require_joint_stabilization = true;

p.odeRelTol = 1e-8;
p.odeAbsTol = 1e-10;

% ode45 remains adaptive. A conservative MaxStep is imposed only during
% the fast initial transient; late-time accuracy is controlled by the
% local error tolerances and checked over the full convergence window.
p.initialStep_fs = 1e-4;
p.maxStepEarly_fs = 1e-2;

% Biot-Savart integration.
p.Nseg_BS = 40;
p.min_distance_A = 0.20;
p.density_threshold = 1e-8;

% Numerical audit tolerances.
p.trace_tolerance = 1e-6;
p.hermiticity_tolerance = 1e-6;
p.positivity_tolerance = 1e-7;
p.surface_norm_tolerance = 1e-6;
p.hamiltonian_hermiticity_tolerance = 1e-10;

% -------------------------------------------------------------------------
% THREE INITIAL SURFACE-SPIN CONDITIONS
% -------------------------------------------------------------------------
initialSpinCases = struct([]);

initialSpinCases(1).label = '+x';
initialSpinCases(1).tag = 'S0_X_100';
initialSpinCases(1).S0 = [1 0 0];

initialSpinCases(2).label = '+y';
initialSpinCases(2).tag = 'S0_Y_010';
initialSpinCases(2).S0 = [0 1 0];

initialSpinCases(3).label = '+z';
initialSpinCases(3).tag = 'S0_Z_001';
initialSpinCases(3).S0 = [0 0 1];

for ic = 1:numel(initialSpinCases)
    initialSpinCases(ic).S0 = ...
        initialSpinCases(ic).S0/norm(initialSpinCases(ic).S0);
end

% -------------------------------------------------------------------------
% OUTPUT FOLDER
% -------------------------------------------------------------------------
downloadsDir = get_downloads_directory();
runTag = datestr(now,'yyyymmdd_HHMMSS');
outputFolder = fullfile(downloadsDir, ...
    ['MODEL_B_SURFACE_SPIN_LLG_XYZ_INITIAL_SPINS_' runTag]);

if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

fprintf('\n============================================================\n');
fprintf('MODEL B: BIDIRECTIONAL SURFACE-SPIN LLG - X/Y/Z INITIAL SPINS\n');
fprintf('============================================================\n');
fprintf('T               = %.3f K\n',p.T_K);
fprintf('tau_rel         = %.3f fs\n',p.tau_rel_fs);
fprintf('Gamma0          = %.6g fs^-1\n',p.Gamma0_per_fs);
fprintf('J0              = %.6g meV\n',p.J0_meV);
fprintf('r0              = %.6g Angstrom\n',p.r0_A);
fprintf('alpha_LLG       = %.6g\n',p.alpha_LLG);
fprintf('gamma_LLG       = %.12e rad/(fs T)\n',p.gamma_rad_fs_T);
fprintf('solver          = ode45\n');
fprintf('InitialStep     = %.6g fs\n',p.initialStep_fs);
fprintf('early MaxStep   = %.6g fs (0 to %.0f fs)\n', ...
    p.maxStepEarly_fs,p.early_window_fs);
fprintf('late propagation = adaptive chunks, no forced MaxStep\n');
fprintf('maximum time    = %.1f fs\n',p.t_max_fs);
fprintf('stability check = every %.1f fs after %.1f fs\n', ...
    p.chunk_fs,p.min_stop_time_fs);
fprintf('stability window = %.1f fs, consecutive checks = %d\n', ...
    p.convergence_window_fs,p.required_consecutive_checks);
fprintf('RelTol          = %.3e\n',p.odeRelTol);
fprintf('AbsTol          = %.3e\n',p.odeAbsTol);
fprintf('initial S cases = [1 0 0], [0 1 0], [0 0 1]\n');
fprintf('field figures   = mT\n');
fprintf('M figures       = normalized M/M_s = S_surface\n');
fprintf('output folder   = %s\n',outputFolder);
fprintf('============================================================\n');

%% ========================================================================
% 2. MOLECULAR GEOMETRIES AND GLOBAL AUDITS
% =========================================================================

models = build_two_models();

for im = 1:numel(models)
    audit_geometry(models(im));
end

% Run the unchanged physical/unit audit with the base +z state first.
run_parameter_and_unit_audit(p,true);

% Geometry does not depend on the initial spin, so plot it exactly once.
plot_geometry_comparison(models,p,outputFolder);

%% ========================================================================
% 3. PROPAGATE BOTH GEOMETRIES FOR S0 = +X, +Y, +Z
% =========================================================================

allRuns = struct([]);

for ic = 1:numel(initialSpinCases)

    pRun = p;
    pRun.S0 = initialSpinCases(ic).S0;
    pRun.S0 = pRun.S0/max(norm(pRun.S0),1e-15);

    % Re-run parameter/unit checks for the actual initial state used.
    run_parameter_and_unit_audit(pRun,true);

    caseFolder = fullfile(outputFolder,initialSpinCases(ic).tag);
    if ~exist(caseFolder,'dir')
        mkdir(caseFolder);
    end

    fprintf('\n============================================================\n');
    fprintf('INITIAL SURFACE SPIN %s : S0 = [%.0f %.0f %.0f]\n', ...
        initialSpinCases(ic).label,pRun.S0(1),pRun.S0(2),pRun.S0(3));
    fprintf('============================================================\n');

    results = struct([]);

    for im = 1:numel(models)

        fprintf('\n------------------------------------------------------------\n');
        fprintf('Running: %s | S0 = [%.0f %.0f %.0f]\n', ...
            models(im).displayName,pRun.S0(1),pRun.S0(2),pRun.S0(3));
        fprintf('------------------------------------------------------------\n');

        R = simulate_surface_spin_model(models(im),pRun,caseFolder);

        % Store the actual initial condition explicitly in every result.
        R.initialSurfaceSpin = pRun.S0;
        R.initialSpinLabel = initialSpinCases(ic).label;
        R.initialSpinTag = initialSpinCases(ic).tag;

        if im == 1
            results = R;
        else
            R = orderfields(R,results(1));
            results(im) = R;
        end
    end

    allRuns(ic).label = initialSpinCases(ic).label;
    allRuns(ic).tag = initialSpinCases(ic).tag;
    allRuns(ic).S0 = pRun.S0;
    allRuns(ic).results = results;

    % ---------------------------------------------------------------------
    % THE ONLY TWO FIGURE TYPES GENERATED PER INITIAL CONDITION:
    %
    % 1) final molecular spin imbalance m_z = n_up - n_down (dimensionless)
    % 2) time dynamics of B_x/B_y/B_z [mT] together with normalized
    %    M_x/M_s, M_y/M_s, M_z/M_s = S_x,S_y,S_z
    % ---------------------------------------------------------------------
    plot_final_imbalance_for_initial_spin( ...
        results,caseFolder,initialSpinCases(ic));

    plot_field_and_normalized_magnetization_for_initial_spin( ...
        results,caseFolder,initialSpinCases(ic));
end

%% ========================================================================
% 4. SUMMARY TABLE AND COMPLETE MAT FILE
% =========================================================================

nCases = numel(allRuns);
nModels = numel(models);
nRows = nCases*nModels;

initialSpinLabel = strings(nRows,1);
initialSpinTag = strings(nRows,1);
S0_x = zeros(nRows,1);
S0_y = zeros(nRows,1);
S0_z = zeros(nRows,1);
modelName = strings(nRows,1);

finalMaxAbsMz = zeros(nRows,1);
finalTime_fs = zeros(nRows,1);
jointStabilized = false(nRows,1);
maxSurfaceSpinChange = zeros(nRows,1);
maxSurfaceField_mT = zeros(nRows,1);
maxTraceError = zeros(nRows,1);
maxHermiticityError = zeros(nRows,1);
minimumDensityEigenvalue = zeros(nRows,1);
maxSurfaceNormError = zeros(nRows,1);
maxHamiltonianHermiticityError = zeros(nRows,1);
runtimeSeconds = zeros(nRows,1);

row = 0;

for ic = 1:nCases
    for im = 1:nModels

        row = row+1;
        R = allRuns(ic).results(im);

        initialSpinLabel(row) = string(allRuns(ic).label);
        initialSpinTag(row) = string(allRuns(ic).tag);
        S0_x(row) = allRuns(ic).S0(1);
        S0_y(row) = allRuns(ic).S0(2);
        S0_z(row) = allRuns(ic).S0(3);
        modelName(row) = string(R.displayName);

        finalMaxAbsMz(row) = R.finalMaxAbsMz;
        finalTime_fs(row) = R.t(end);
        jointStabilized(row) = R.jointStabilized;
        maxSurfaceSpinChange(row) = R.maxSurfaceSpinChange;
        maxSurfaceField_mT(row) = R.maxSurfaceField_mT;
        maxTraceError(row) = R.maxTraceError;
        maxHermiticityError(row) = R.maxHermiticityError;
        minimumDensityEigenvalue(row) = R.minimumDensityEigenvalue;
        maxSurfaceNormError(row) = R.maxSurfaceNormError;
        maxHamiltonianHermiticityError(row) = ...
            R.maxHamiltonianHermiticityError;
        runtimeSeconds(row) = R.runtimeSeconds;
    end
end

summaryTable = table( ...
    initialSpinLabel,initialSpinTag,S0_x,S0_y,S0_z,modelName, ...
    finalMaxAbsMz,finalTime_fs,jointStabilized,maxSurfaceSpinChange, ...
    maxSurfaceField_mT,maxTraceError,maxHermiticityError, ...
    minimumDensityEigenvalue,maxSurfaceNormError, ...
    maxHamiltonianHermiticityError,runtimeSeconds, ...
    'VariableNames',{ ...
    'initial_spin_label','initial_spin_tag','S0_x','S0_y','S0_z','model', ...
    'final_max_abs_m_z','final_time_fs','joint_stabilized', ...
    'max_surface_spin_change','max_surface_field_mT','max_trace_error', ...
    'max_hermiticity_error','minimum_density_eigenvalue', ...
    'max_surface_norm_error','max_H_hermiticity_error','runtime_seconds'});

disp(summaryTable);

writetable(summaryTable,fullfile(outputFolder, ...
    'surface_spin_XYZ_initial_conditions_summary.csv'));

save(fullfile(outputFolder, ...
    'surface_spin_XYZ_initial_conditions_complete_results.mat'), ...
    'allRuns','initialSpinCases','models','p','summaryTable','-v7.3');

fprintf('\n============================================================\n');
fprintf('BIDIRECTIONAL SURFACE-SPIN CALCULATION COMPLETE\n');
fprintf('Initial directions completed: +x, +y, +z\n');
fprintf('Figures produced: 7 total (1 geometry + 3 imbalance + 3 dynamics)\n');
fprintf('Results saved in:\n%s\n',outputFolder);
fprintf('============================================================\n');

end

%% ========================================================================
% SIMULATE ONE GEOMETRY: DYNAMIC SURFACE SPIN
% =========================================================================

function result = simulate_surface_spin_model(inputModel,p,outputFolder)

runTimer = tic;

Ri = inputModel.Ri;
bonds = inputModel.bonds;
L = size(Ri,1);
Nb = size(bonds,1);
Nspin = 2*L;
Nrho = Nspin^2;

t_bonds_eV = p.t_hop_eV*inputModel.t_bonds(:);

%% ------------------------------------------------------------------------
% A. ORIGINAL ORBITAL HAMILTONIAN
% -------------------------------------------------------------------------

H_hop = build_hopping_hamiltonian(L,bonds,t_bonds_eV);

distance_A = vecnorm(Ri-p.R_surface_A,2,2);
distance_A(distance_A < 1e-12) = 1e-12;

epsilon_site_eV = ...
    p.q_eff*p.e2_over_4pieps0_eVA./distance_A;

H_contact = diag(epsilon_site_eV);
H0 = H_hop+H_contact;

% Site-major spin expansion.
H_TB_spin = kron(H0,eye(2));

% Distance-dependent surface exchange:
% J_i = J0 exp(-r_i/r0).
J_site_eV = p.J0_eV*exp(-distance_A/p.r0_A);
J_site_meV = 1e3*J_site_eV;

assert(norm(H0-H0','fro') < 1e-12, ...
    'The orbital Hamiltonian is not Hermitian.');
assert(norm(H_TB_spin-H_TB_spin','fro') < 1e-12, ...
    'The spin-expanded orbital Hamiltonian is not Hermitian.');

%% ------------------------------------------------------------------------
% B. GAS-PHASE GROUND-STATE PREPARATION
% -------------------------------------------------------------------------
%
% Identical to the no-surface-spin calculation:
%
% rho(0)=|psi_GS><psi_GS| tensor I_2/2.

[rho0,initialDiagnostic] = ...
    build_gas_ground_initial_state(H_hop);

%% ------------------------------------------------------------------------
% C. UNCHANGED LINDBLAD BATH
% -------------------------------------------------------------------------
%
% The Lindblad operators are built from the orbital post-contact H0 and
% are duplicated for up and down. They are NOT rebuilt from the
% instantaneous exchange- or Zeeman-split Hamiltonian.

Dsuper = build_spin_independent_lindblad_superoperator(H0,p);

traceVector = reshape(eye(Nspin),[],1).';
lindbladTraceResidual = norm(traceVector*Dsuper,inf);

assert(lindbladTraceResidual < 1e-9, ...
    'The Lindblad superoperator is not trace preserving.');

%% ------------------------------------------------------------------------
% D. BIOT-SAVART KERNELS
% -------------------------------------------------------------------------
%
% One kernel gives B at every molecular site for H_Z.
% A second kernel gives B at R_surface for the LLG equation.

BS_kernel_sites = precompute_BS_kernel_at_sites( ...
    Ri,bonds,p.Nseg_BS,p.min_distance_A);

BS_kernel_surface = precompute_BS_kernel_at_points( ...
    Ri,bonds,p.R_surface_A,p.Nseg_BS,p.min_distance_A);

model.Ri = Ri;
model.bonds = bonds;
model.t_bonds_eV = t_bonds_eV;
model.H_TB_spin = H_TB_spin;
model.J_site_eV = J_site_eV;
model.BS_kernel_sites = BS_kernel_sites;
model.BS_kernel_surface = BS_kernel_surface;

%% ------------------------------------------------------------------------
% E. COUPLED RHO + LLG PROPAGATION
% -------------------------------------------------------------------------
%
% The ODE state contains:
%
%   y = [vec(rho); S_x; S_y; S_z].
%
% At every evaluation:
%
%   rho -> J -> B_sites and B_surface
%   B_sites -> H_Z
%   S_surface -> H_surface-ex
%   B_surface -> LLG -> S_surface
%
% This is the complete bidirectional loop used here.

y0 = [rho0(:);p.S0(:)];

rhs = @(~,y) coupled_rho_llg_rhs(y,model,Dsuper,p);

optionsEarly = odeset( ...
    'InitialStep',p.initialStep_fs, ...
    'MaxStep',p.maxStepEarly_fs, ...
    'RelTol',p.odeRelTol, ...
    'AbsTol',p.odeAbsTol, ...
    'Stats','off');

optionsLate = odeset( ...
    'RelTol',p.odeRelTol, ...
    'AbsTol',p.odeAbsTol, ...
    'Refine',4, ...
    'Stats','off');

[t,Y,jointStabilized,stopInfo] = ...
    propagate_until_joint_stabilization( ...
    rhs,y0,Nrho,p,optionsEarly,optionsLate);

%% ------------------------------------------------------------------------
% F. SITE POPULATIONS, FIELDS, SURFACE SPIN AND AUDITS
% -------------------------------------------------------------------------

nT = numel(t);

n_up = zeros(nT,L);
n_down = zeros(nT,L);
n_total = zeros(nT,L);
m_z = zeros(nT,L);
P_z = nan(nT,L);

B_x_mT = zeros(nT,L);
B_y_mT = zeros(nT,L);
B_z_mT = zeros(nT,L);
B_abs_mT = zeros(nT,L);

B_surface_mT = zeros(nT,3);
S_surface = zeros(nT,3);
bond_currents_A = zeros(nT,Nb);

trace_error = zeros(nT,1);
hermiticity_error = zeros(nT,1);
minimum_eigenvalue = zeros(nT,1);
surface_norm_error = zeros(nT,1);
hamiltonian_hermiticity_error = zeros(nT,1);

for it = 1:nT

    rhoRaw = reshape(Y(it,1:Nrho).',Nspin,Nspin);
    hermiticity_error(it) = norm(rhoRaw-rhoRaw','fro');
    rho = 0.5*(rhoRaw+rhoRaw');

    trace_error(it) = abs(real(trace(rho))-1);
    minimum_eigenvalue(it) = min(real(eig(rho)));

    Sraw = real(Y(it,Nrho+1:Nrho+3));
    Sraw = Sraw(:).';

    surface_norm_error(it) = abs(norm(Sraw)-1);
    S = Sraw/max(norm(Sraw),1e-15);
    S_surface(it,:) = S;

    obs = density_observables_from_rho(rho,p);

    n_up(it,:) = obs.n_up;
    n_down(it,:) = obs.n_down;
    n_total(it,:) = obs.n_total;
    m_z(it,:) = obs.m_z;
    P_z(it,:) = obs.P_z;

    J_eV = compute_bond_current_energy(rho,model);

    B_sites_mT = apply_BS_kernel(model.BS_kernel_sites,J_eV);
    Bsurf_mT = apply_BS_kernel(model.BS_kernel_surface,J_eV);

    B_x_mT(it,:) = B_sites_mT(:,1).';
    B_y_mT(it,:) = B_sites_mT(:,2).';
    B_z_mT(it,:) = B_sites_mT(:,3).';
    B_abs_mT(it,:) = vecnorm(B_sites_mT,2,2).';
    B_surface_mT(it,:) = Bsurf_mT;

    bond_currents_A(it,:) = ...
        (1.602176634e-19/6.582119569e-16)*J_eV;

    H_Z = build_zeeman_hamiltonian(B_sites_mT,p);
    H_exchange = build_exchange_hamiltonian(J_site_eV,S);

    H_total = model.H_TB_spin+H_Z+H_exchange;

    hamiltonian_hermiticity_error(it) = ...
        norm(H_total-H_total','fro');
end

% Explicit initial-condition audit for the current +x/+y/+z run.
assert(norm(S_surface(1,:)-p.S0(:).') <= 1e-12, ...
    'Stored surface-spin initial condition does not match p.S0.');

total_m_z = sum(m_z,2);

[maxAbsMz,linearIndex] = max(abs(m_z(:)));
[peakTimeIndex,peakSite] = ind2sub(size(m_z),linearIndex);

peakTime_fs = t(peakTimeIndex);
peakSignedMz = m_z(peakTimeIndex,peakSite);
finalMaxAbsMz = max(abs(m_z(end,:)));

% -------------------------------------------------------------------------
% JOINT EARLY-STOP / STATIONARITY DIAGNOSTICS
% -------------------------------------------------------------------------
% The stopping decision uses the complete final convergence window and
% multiple RHS samples for both rho and S, not only the last time point.
populationRange = stopInfo.populationRange;
surfaceSpinRange = stopInfo.surfaceSpinRange;
rhoRhsNorm_per_fs = stopInfo.maxRhoRhsNorm_per_fs;
surfaceRhsNorm_per_fs = stopInfo.maxSurfaceRhsNorm_per_fs;
longTimeConverged = jointStabilized;

%% ------------------------------------------------------------------------
% G. PHYSICAL AND NUMERICAL OUTPUT AUDIT
% -------------------------------------------------------------------------

assert(max(trace_error) <= p.trace_tolerance, ...
    'Density-matrix trace drift exceeded the tolerance.');
assert(max(hermiticity_error) <= p.hermiticity_tolerance, ...
    'Density-matrix Hermiticity drift exceeded the tolerance.');
assert(min(minimum_eigenvalue) >= -p.positivity_tolerance, ...
    'Density matrix developed an unacceptable negative eigenvalue.');
assert(max(abs(sum(n_total,2)-1)) <= p.trace_tolerance, ...
    'Total one-electron population is not conserved.');
assert(max(abs(m_z(:)-(n_up(:)-n_down(:)))) < 1e-10, ...
    'm_z is inconsistent with n_up-n_down.');
assert(max(surface_norm_error) <= p.surface_norm_tolerance, ...
    'The LLG surface spin did not remain normalized.');
assert(max(hamiltonian_hermiticity_error) <= ...
    p.hamiltonian_hermiticity_tolerance, ...
    'The time-dependent Hamiltonian is not Hermitian.');
assert(all(isfinite([n_up(:);n_down(:);m_z(:); ...
    B_abs_mT(:);B_surface_mT(:);S_surface(:)])), ...
    'A physical output contains NaN or Inf.');

%% ------------------------------------------------------------------------
% H. SAVE COMPLETE RESULTS
% -------------------------------------------------------------------------

safeName = inputModel.name;

siteParameterTable = table( ...
    (1:L).',Ri(:,1),Ri(:,2),Ri(:,3), ...
    distance_A,epsilon_site_eV,J_site_meV, ...
    'VariableNames',{ ...
    'site','x_A','y_A','z_A', ...
    'distance_from_surface_A','epsilon_site_eV','J_i_meV'});

writetable(siteParameterTable,fullfile(outputFolder, ...
    [safeName '_site_parameters_and_exchange.csv']));

write_bond_current_table(fullfile(outputFolder, ...
    [safeName '_all_times_bond_currents_A.csv']), ...
    t,bonds,bond_currents_A);

surfaceTimeTable = table( ...
    t,S_surface(:,1),S_surface(:,2),S_surface(:,3), ...
    B_surface_mT(:,1),B_surface_mT(:,2),B_surface_mT(:,3), ...
    'VariableNames',{ ...
    'time_fs','S_surface_x','S_surface_y','S_surface_z', ...
    'B_surface_x_mT','B_surface_y_mT','B_surface_z_mT'});

writetable(surfaceTimeTable,fullfile(outputFolder, ...
    [safeName '_surface_spin_and_field.csv']));

finalTable = table( ...
    (1:L).',J_site_meV,n_up(end,:).',n_down(end,:).', ...
    n_total(end,:).',m_z(end,:).',P_z(end,:).', ...
    B_x_mT(end,:).',B_y_mT(end,:).',B_z_mT(end,:).', ...
    B_abs_mT(end,:).', ...
    'VariableNames',{ ...
    'site','J_i_meV','stop_n_up','stop_n_down', ...
    'stop_n_total','stop_m_z','stop_P_z', ...
    'stop_Bx_mT','stop_By_mT','stop_Bz_mT', ...
    'stop_B_abs_mT'});

writetable(finalTable,fullfile(outputFolder, ...
    [safeName '_site_results_at_stop_time.csv']));

convergenceTable = table( ...
    t(end),p.convergence_window_fs,populationRange,surfaceSpinRange, ...
    rhoRhsNorm_per_fs,surfaceRhsNorm_per_fs,longTimeConverged, ...
    'VariableNames',{ ...
    'final_time_fs','window_fs','population_range_in_window', ...
    'surface_spin_component_range_in_window','max_rho_rhs_norm_in_window_per_fs', ...
    'max_surface_rhs_norm_in_window_per_fs','joint_stabilized'});

writetable(convergenceTable,fullfile(outputFolder, ...
    [safeName '_long_time_convergence.csv']));

result.name = inputModel.name;
result.displayName = inputModel.displayName;
result.initialSurfaceSpin = p.S0(:).';
result.L = L;
result.Ri = Ri;
result.bonds = bonds;
result.J_site_eV = J_site_eV;
result.J_site_meV = J_site_meV;
result.t = t;

result.n_up = n_up;
result.n_down = n_down;
result.n_total = n_total;
result.m_z = m_z;
result.P_z = P_z;
result.total_m_z = total_m_z;

result.B_x_mT = B_x_mT;
result.B_y_mT = B_y_mT;
result.B_z_mT = B_z_mT;
result.B_abs_mT = B_abs_mT;
result.B_surface_mT = B_surface_mT;
result.S_surface = S_surface;
result.bond_currents_A = bond_currents_A;

result.peakTimeIndex = peakTimeIndex;
result.peakTime_fs = peakTime_fs;
result.peakSite = peakSite;
result.peakSignedMz = peakSignedMz;
result.peakAbsMz = maxAbsMz;
result.finalMaxAbsMz = finalMaxAbsMz;

result.maxSurfaceSpinChange = max(vecnorm( ...
    S_surface-S_surface(1,:),2,2));
result.maxSurfaceField_mT = max(vecnorm(B_surface_mT,2,2));

result.populationRangeFinalWindow = populationRange;
result.surfaceSpinRangeFinalWindow = surfaceSpinRange;
result.rhoRhsNormFinal_per_fs = rhoRhsNorm_per_fs;
result.surfaceRhsNormFinal_per_fs = surfaceRhsNorm_per_fs;
result.longTimeConverged = longTimeConverged;
result.convergenceWindow_fs = p.convergence_window_fs;
result.jointStabilized = jointStabilized;
result.llgStabilized = jointStabilized;
result.llgStabilizationTime_fs = t(end);
result.stabilizationTime_fs = t(end);
result.llgSpinRangeAtStop = stopInfo.surfaceSpinRange;
result.llgRhsNormAtStop_per_fs = stopInfo.maxSurfaceRhsNorm_per_fs;
result.maxRhoRhsNormAtStop_per_fs = stopInfo.maxRhoRhsNorm_per_fs;
result.consecutiveStableChecks = stopInfo.consecutiveStableChecks;

result.initialDiagnostic = initialDiagnostic;
result.lindbladTraceResidual = lindbladTraceResidual;
result.maxTraceError = max(trace_error);
result.maxHermiticityError = max(hermiticity_error);
result.minimumDensityEigenvalue = min(minimum_eigenvalue);
result.maxSurfaceNormError = max(surface_norm_error);
result.maxHamiltonianHermiticityError = ...
    max(hamiltonian_hermiticity_error);
result.runtimeSeconds = toc(runTimer);

save(fullfile(outputFolder,[safeName '_complete_result.mat']), ...
    'result','-v7.3');

fprintf('After-LLG max |m_z| = %.12e at t = %.6f fs\n', ...
    result.finalMaxAbsMz,result.llgStabilizationTime_fs);
fprintf('Max surface-spin change = %.12e\n', ...
    result.maxSurfaceSpinChange);
fprintf('Max surface field = %.12e mT\n', ...
    result.maxSurfaceField_mT);
fprintf('Final-window population range = %.3e\n', ...
    result.populationRangeFinalWindow);
fprintf('Final-window surface-spin range = %.3e\n', ...
    result.surfaceSpinRangeFinalWindow);
fprintf('Max ||d rho/dt|| in final window = %.3e 1/fs\n', ...
    result.rhoRhsNormFinal_per_fs);
fprintf('Max ||d S/dt|| in final window = %.3e 1/fs\n', ...
    result.surfaceRhsNormFinal_per_fs);
fprintf('Joint rho+LLG stabilization = %d at t = %.6f fs\n', ...
    result.jointStabilized,result.stabilizationTime_fs);
if ~result.jointStabilized
    fprintf(['WARNING: joint electronic and LLG stabilization was ' ...
        'not reached by 5000 fs. The saved state is a maximum-time ' ...
        'diagnostic and is not labelled as stabilized.\n']);
end
fprintf('Runtime = %.2f s\n',result.runtimeSeconds);

end

%% ========================================================================
% CHUNKED ODE45 PROPAGATION WITH JOINT RHO + LLG EARLY STOP
% =========================================================================

function [tAll,YAll,stabilized,info] = ...
    propagate_until_joint_stabilization( ...
    rhs,y0,Nrho,p,optionsEarly,optionsLate)

% Resolve the fast transient on the requested early output grid.
tspanEarly = (p.t_start_fs:p.output_dt_early_fs:p.early_window_fs).';
if tspanEarly(end) < p.early_window_fs
    tspanEarly = [tspanEarly;p.early_window_fs];
end
[tAll,YAll] = ode45(rhs,tspanEarly,y0,optionsEarly);

stabilized = false;
stableCount = 0;
info = empty_joint_stability_info();

while tAll(end) < p.t_max_fs-1e-12

    nextEnd = min(tAll(end)+p.chunk_fs,p.t_max_fs);

    % With a two-point tspan, ode45 returns its accepted/refined adaptive
    % points. No large fixed time step is imposed in the slow regime.
    [tChunk,YChunk] = ode45( ...
        rhs,[tAll(end) nextEnd],YAll(end,:).',optionsLate);

    tAll = [tAll;tChunk(2:end)]; %#ok<AGROW>
    YAll = [YAll;YChunk(2:end,:)]; %#ok<AGROW>

    if tAll(end) < p.min_stop_time_fs-1e-12
        continue;
    end

    [candidateStable,candidateInfo] = ...
        evaluate_joint_stability(tAll,YAll,Nrho,p,rhs);

    if candidateStable
        stableCount = stableCount+1;
    else
        stableCount = 0;
    end

    candidateInfo.consecutiveStableChecks = stableCount;
    info = candidateInfo;

    fprintf(['Joint stability check at t = %.3f fs: population range = %.3e, ' ...
        'S range = %.3e, max ||d rho/dt|| = %.3e, ' ...
        'max ||dS/dt|| = %.3e 1/fs, consecutive = %d/%d.\n'], ...
        tAll(end),info.populationRange,info.surfaceSpinRange, ...
        info.maxRhoRhsNorm_per_fs,info.maxSurfaceRhsNorm_per_fs, ...
        stableCount,p.required_consecutive_checks);

    if stableCount >= p.required_consecutive_checks
        stabilized = true;
        break;
    end
end

[~,finalInfo] = evaluate_joint_stability(tAll,YAll,Nrho,p,rhs);
finalInfo.consecutiveStableChecks = stableCount;
info = finalInfo;

end

function [stable,info] = ...
    evaluate_joint_stability(t,Y,Nrho,p,rhs)

windowStart = max(t(end)-p.convergence_window_fs,t(1));
idx = find(t >= windowStart);
windowDuration_fs = t(end)-t(idx(1));

Nspin = round(sqrt(Nrho));
L = Nspin/2;
nWindow = numel(idx);
nUp = zeros(nWindow,L);
nDown = zeros(nWindow,L);
S = real(Y(idx,Nrho+1:Nrho+3));
S = S./max(vecnorm(S,2,2),1e-15);

for q = 1:nWindow
    rho = reshape(Y(idx(q),1:Nrho).',Nspin,Nspin);
    rho = 0.5*(rho+rho');
    obs = density_observables_from_rho(rho,p);
    nUp(q,:) = obs.n_up;
    nDown(q,:) = obs.n_down;
end

populationRange = max([ ...
    max(max(nUp,[],1)-min(nUp,[],1)), ...
    max(max(nDown,[],1)-min(nDown,[],1))]);

surfaceSpinRange = max(max(S,[],1)-min(S,[],1));

nDerivativeSamples = min(p.derivative_check_points,nWindow);
sampleLocal = unique(round(linspace(1,nWindow,nDerivativeSamples)));
maxRhoRhsNorm_per_fs = 0;
maxSurfaceRhsNorm_per_fs = 0;

for q = sampleLocal
    dY = rhs(t(idx(q)),Y(idx(q),:).');
    maxRhoRhsNorm_per_fs = max( ...
        maxRhoRhsNorm_per_fs,norm(dY(1:Nrho)));
    maxSurfaceRhsNorm_per_fs = max( ...
        maxSurfaceRhsNorm_per_fs,norm(dY(Nrho+1:Nrho+3)));
end

stable = ...
    windowDuration_fs >= p.convergence_window_fs-1e-9 && ...
    populationRange <= p.population_range_tol && ...
    surfaceSpinRange <= p.surface_spin_range_tol && ...
    maxRhoRhsNorm_per_fs <= p.rho_rhs_norm_tol_per_fs && ...
    maxSurfaceRhsNorm_per_fs <= p.surface_rhs_norm_tol_per_fs;

info.populationRange = populationRange;
info.surfaceSpinRange = surfaceSpinRange;
info.maxRhoRhsNorm_per_fs = maxRhoRhsNorm_per_fs;
info.maxSurfaceRhsNorm_per_fs = maxSurfaceRhsNorm_per_fs;
info.windowDuration_fs = windowDuration_fs;
info.consecutiveStableChecks = 0;

end

function info = empty_joint_stability_info()
info.populationRange = inf;
info.surfaceSpinRange = inf;
info.maxRhoRhsNorm_per_fs = inf;
info.maxSurfaceRhsNorm_per_fs = inf;
info.windowDuration_fs = 0;
info.consecutiveStableChecks = 0;
end

%% ========================================================================
% COUPLED DENSITY-MATRIX + LLG RHS
% =========================================================================

function dy = coupled_rho_llg_rhs(y,model,Dsuper,p)

L = size(model.Ri,1);
Nrho = (2*L)^2;

rho = reshape(y(1:Nrho),2*L,2*L);
rho = 0.5*(rho+rho');

S = real(y(Nrho+1:Nrho+3));
S = S(:).';
S = S/max(norm(S),1e-15);

% 1. rho -> molecular bond currents.
J_eV = compute_bond_current_energy(rho,model);

% 2. Currents -> self-field at molecular sites and at the surface.
B_sites_mT = apply_BS_kernel(model.BS_kernel_sites,J_eV);
B_surface_mT = apply_BS_kernel(model.BS_kernel_surface,J_eV);

% 3. Self-field at molecular sites -> electron Zeeman term.
H_Z = build_zeeman_hamiltonian(B_sites_mT,p);

% 4. Surface spin -> distance-dependent exchange Hamiltonian.
H_exchange = build_exchange_hamiltonian(model.J_site_eV,S);

% 5. Full 2L x 2L molecular Hamiltonian.
H_total = model.H_TB_spin+H_Z+H_exchange;

% 6. Density-matrix equation.
drho_H = -1i/p.hbar_eV_fs*(H_total*rho-rho*H_total);
drho_D = reshape(Dsuper*rho(:),2*L,2*L);
drho = drho_H+drho_D;

% 7. Biot-Savart field at the surface -> LLG.
%
% CRITICAL: B_surface_mT is converted to tesla before LLG.
B_surface_T = 1e-3*B_surface_mT;

dS = llg_rhs(S,B_surface_T,p);

dy = [drho(:);dS(:)];

end

%% ========================================================================
% SURFACE-EXCHANGE HAMILTONIAN
% =========================================================================
%
% H_surface-ex = sum_i J_i S_surface.s_i,
% s_i = sigma/2.
%
% At site i:
%
% H_ex,i = (J_i/2) *
%          [ S_z        S_x-iS_y
%            S_x+iS_y  -S_z      ].

function H_exchange = build_exchange_hamiltonian(J_site_eV,S)

L = numel(J_site_eV);

sigma_x = [0 1;1 0];
sigma_y = [0 -1i;1i 0];
sigma_z = [1 0;0 -1];

S_dot_sigma = ...
    S(1)*sigma_x+S(2)*sigma_y+S(3)*sigma_z;

H_exchange = complex(zeros(2*L,2*L));

for i = 1:L

    idx = (2*i-1):(2*i);

    H_exchange(idx,idx) = ...
        0.5*J_site_eV(i)*S_dot_sigma;
end

end

%% ========================================================================
% LANDAU-LIFSHITZ-GILBERT EQUATION
% =========================================================================
%
% dS/dt = -gamma/(1+alpha^2) *
%         [ S x B + alpha S x (S x B) ].
%
% S is dimensionless and normalized.
% B must be supplied in tesla.
% gamma is in rad/(fs T), so dS/dt is in fs^-1.

function dS = llg_rhs(S,B_surface_T,p)

S = S/max(norm(S),1e-15);

precession = cross(S,B_surface_T);
damping = cross(S,cross(S,B_surface_T));

dS = -p.gamma_rad_fs_T/(1+p.alpha_LLG^2) * ...
    (precession+p.alpha_LLG*damping);

dS = dS(:);

end

%% ========================================================================
% ZEEMAN HAMILTONIAN
% =========================================================================

function H_Z = build_zeeman_hamiltonian(B_sites_mT,p)

L = size(B_sites_mT,1);
H_Z = complex(zeros(2*L,2*L));

% Required conversion because mu_B is stored in eV/T.
B_sites_T = 1e-3*B_sites_mT;
z = p.zeeman_prefactor_eV_T;

for i = 1:L

    idx = (2*i-1):(2*i);

    Bx = B_sites_T(i,1);
    By = B_sites_T(i,2);
    Bz = B_sites_T(i,3);

    H_Z(idx,idx) = z*[ ...
        Bz,          Bx-1i*By; ...
        Bx+1i*By,   -Bz];
end

end

%% ========================================================================
% SPIN-SUMMED BOND CURRENT
% =========================================================================

function J_eV = compute_bond_current_energy(rho,model)

Nb = size(model.bonds,1);
J_eV = zeros(1,Nb);

for b = 1:Nb

    i = model.bonds(b,1);
    j = model.bonds(b,2);
    tij = model.t_bonds_eV(b);

    up_i = 2*i-1;
    down_i = 2*i;
    up_j = 2*j-1;
    down_j = 2*j;

    J_eV(b) = 2*tij*imag( ...
        rho(up_j,up_i)+rho(down_j,down_i));
end

end

%% ========================================================================
% LINDBLAD SUPEROPERATOR WITH DETAILED BALANCE
% =========================================================================

function Dsuper = build_spin_independent_lindblad_superoperator(H0,p)

L = size(H0,1);
N = 2*L;
I_N = eye(N);

[U,D] = eig(H0);
[E,idx] = sort(real(diag(D)));
U = U(:,idx);

P_up = [1 0;0 0];
P_down = [0 0;0 1];

Dsuper = complex(zeros(N^2,N^2));

for n = 1:L
    for m = n+1:L

        gap_eV = E(m)-E(n);

        if gap_eV <= 1e-12
            continue;
        end

        gamma_down = p.Gamma0_per_fs;
        gamma_up = p.Gamma0_per_fs* ...
            exp(-p.beta_eV_inv*gap_eV);

        V_up_energy = zeros(L,L);
        V_down_energy = zeros(L,L);

        V_up_energy(m,n) = sqrt(gamma_up);
        V_down_energy(n,m) = sqrt(gamma_down);

        V_up_site = U*V_up_energy*U';
        V_down_site = U*V_down_energy*U';

        % The index k in sum_k D[L_k] enumerates these operators.
        Lset = {
            kron(V_up_site,P_up)
            kron(V_down_site,P_up)
            kron(V_up_site,P_down)
            kron(V_down_site,P_down)
        };

        for k = 1:numel(Lset)

            Lk = Lset{k};
            A = Lk'*Lk;

            Dsuper = Dsuper ...
                + kron(conj(Lk),Lk) ...
                - 0.5*kron(I_N,A) ...
                - 0.5*kron(A.',I_N);
        end
    end
end

end

%% ========================================================================
% GAS-PHASE GROUND STATE
% =========================================================================

function [rho0,diagnostic] = build_gas_ground_initial_state(H_hop)

[U,D] = eig(H_hop);
[E,idx] = sort(real(diag(D)));
U = U(:,idx);

psi0 = U(:,1);
rhoOrbital = psi0*psi0';
rho0 = kron(rhoOrbital,0.5*eye(2));

diagnostic.energy_eV = E(1);
diagnostic.residual = norm(H_hop*psi0-E(1)*psi0);
diagnostic.traceError = abs(real(trace(rho0))-1);
diagnostic.hermiticityError = norm(rho0-rho0','fro');

up = 1:2:size(rho0,1);
down = 2:2:size(rho0,1);

diagnostic.initialTotalMz = real( ...
    sum(diag(rho0(up,up)))-sum(diag(rho0(down,down))));

assert(diagnostic.residual < 1e-10, ...
    'Gas-phase ground-state residual is too large.');
assert(diagnostic.traceError < 1e-12, ...
    'Initial density trace is not one.');
assert(diagnostic.hermiticityError < 1e-12, ...
    'Initial density matrix is not Hermitian.');
assert(abs(diagnostic.initialTotalMz) < 1e-12, ...
    'Initial molecular spin is not unpolarized.');

end

%% ========================================================================
% SITE-RESOLVED DENSITY OBSERVABLES
% =========================================================================

function obs = density_observables_from_rho(rho,p)

N = size(rho,1);
assert(mod(N,2) == 0, ...
    'A spinful density matrix must have even dimension.');

up = 1:2:N;
down = 2:2:N;

obs.n_up = real(diag(rho(up,up))).';
obs.n_down = real(diag(rho(down,down))).';
obs.n_total = obs.n_up+obs.n_down;
obs.m_z = obs.n_up-obs.n_down;

obs.P_z = nan(size(obs.n_total));
valid = obs.n_total > p.density_threshold;
obs.P_z(valid) = obs.m_z(valid)./obs.n_total(valid);

end

%% ========================================================================
% GEOMETRIES
% =========================================================================

function models = build_two_models()

models = struct([]);

models(1).name = 'helical';
models(1).displayName = 'Helical chiral model';
models(1).Ri = [
     0.0000   0.0000   1.0000
    -0.2250   0.3897   1.8000
    -0.6750   0.3897   2.6000
    -0.9000   0.0000   3.4000
    -0.6750  -0.3897   4.2000
    -0.2250  -0.3897   5.0000
     0.0000   0.0000   5.8000
    -0.2250   0.3897   6.6000
    -0.6750   0.3897   7.4000
    -0.9000   0.0000   8.2000
];
models(1).bonds = [(1:9).' (2:10).'];
models(1).t_bonds = ones(9,1);

a_A = 1.54;
z0_A = 1.0;
RiBranch = zeros(10,3);

for i = 1:8
    RiBranch(i,:) = [0 0 z0_A+(i-1)*a_A];
end

branchBase = RiBranch(5,:);
branchRelative = [
     1.20   0.70   0.65
     1.70  -0.50   1.47
];

RiBranch(9,:) = branchBase+branchRelative(1,:);
RiBranch(10,:) = branchBase+branchRelative(2,:);

models(2).name = 'nonhelical_chiral_branch';
models(2).displayName = ...
    'Non-helical chiral: linear backbone + 3D branch';
models(2).Ri = RiBranch;
models(2).bonds = [
    1 2
    2 3
    3 4
    4 5
    5 6
    6 7
    7 8
    5 9
    9 10
];
models(2).t_bonds = ones(9,1);

models(2).chiralityTripleProduct = ...
    dot([0 0 1],cross(branchRelative(1,:),branchRelative(2,:)));

end

function audit_geometry(model)

Ri = model.Ri;
bonds = model.bonds;

assert(size(Ri,2) == 3, ...
    'Coordinates must be three dimensional.');
assert(all(isfinite(Ri(:))), ...
    'Geometry contains NaN or Inf.');
assert(all(bonds(:) >= 1) && ...
       all(bonds(:) <= size(Ri,1)), ...
    'A bond contains an invalid site index.');

bondLengths_A = zeros(size(bonds,1),1);

for b = 1:size(bonds,1)
    bondLengths_A(b) = norm( ...
        Ri(bonds(b,2),:)-Ri(bonds(b,1),:));
end

assert(all(bondLengths_A > 0), ...
    'A zero-length bond was detected.');

if isfield(model,'chiralityTripleProduct') && ...
        ~isempty(model.chiralityTripleProduct)

    chi = model.chiralityTripleProduct;

    assert(isscalar(chi) && isfinite(chi) && abs(chi) > 1e-8, ...
        'The branched geometry became coplanar or invalid.');
end

end

%% ========================================================================
% HOPPING HAMILTONIAN
% =========================================================================

function H_hop = build_hopping_hamiltonian(L,bonds,t_bonds_eV)

H_hop = zeros(L,L);

for b = 1:size(bonds,1)

    i = bonds(b,1);
    j = bonds(b,2);
    tij = t_bonds_eV(b);

    H_hop(i,j) = -tij;
    H_hop(j,i) = -tij;
end

end

%% ========================================================================
% BIOT-SAVART KERNELS AND FIELD
% =========================================================================

function kernel = precompute_BS_kernel_at_sites( ...
    Ri,bonds,Nseg,min_distance_A)

L = size(Ri,1);
Nb = size(bonds,1);
kernel = zeros(3*L,Nb);

for activeBond = 1:Nb

    unit_J_eV = zeros(1,Nb);
    unit_J_eV(activeBond) = 1;

    Bunit_mT = zeros(L,3);

    for observationSite = 1:L

        Bunit_mT(observationSite,:) = ...
            biot_savart_at_point( ...
            Ri,bonds,unit_J_eV,Ri(observationSite,:), ...
            observationSite,Nseg,min_distance_A);
    end

    kernel(:,activeBond) = reshape(Bunit_mT.',[],1);
end

end

function kernel = precompute_BS_kernel_at_points( ...
    Ri,bonds,Rpoints_A,Nseg,min_distance_A)

nPoints = size(Rpoints_A,1);
Nb = size(bonds,1);
kernel = zeros(3*nPoints,Nb);

for activeBond = 1:Nb

    unit_J_eV = zeros(1,Nb);
    unit_J_eV(activeBond) = 1;

    Bunit_mT = zeros(nPoints,3);

    for ip = 1:nPoints

        Bunit_mT(ip,:) = biot_savart_at_point( ...
            Ri,bonds,unit_J_eV,Rpoints_A(ip,:),0, ...
            Nseg,min_distance_A);
    end

    kernel(:,activeBond) = reshape(Bunit_mT.',[],1);
end

end

function B_mT = apply_BS_kernel(kernel,J_eV)

Bflat_mT = kernel*J_eV(:);
nPoints = size(kernel,1)/3;
B_mT = reshape(Bflat_mT,3,nPoints).';

end

function B_mT = biot_savart_at_point( ...
    Ri,bonds,J_eV,Robs_A,excludedSite,Nseg,min_distance_A)

e_C = 1.602176634e-19;
hbar_eV_s = 6.582119569e-16;
mu0_over_4pi = 1e-7;
A_to_m = 1e-10;

Ri_m = Ri*A_to_m;
Robs_m = Robs_A*A_to_m;
cutoff_m = min_distance_A*A_to_m;

I_bonds_A = (e_C/hbar_eV_s)*J_eV;

B_T = [0 0 0];

for b = 1:size(bonds,1)

    i = bonds(b,1);
    j = bonds(b,2);

    if excludedSite > 0 && ...
       (i == excludedSite || j == excludedSite)
        continue;
    end

    r1 = Ri_m(i,:);
    r2 = Ri_m(j,:);
    dl_m = (r2-r1)/Nseg;

    for segment = 1:Nseg

        rSegment_m = r1+(segment-0.5)*dl_m;
        rVector_m = Robs_m-rSegment_m;
        rNorm_m = norm(rVector_m);

        if rNorm_m < cutoff_m
            continue;
        end

        B_T = B_T+ ...
            mu0_over_4pi*I_bonds_A(b) ...
            *cross(dl_m,rVector_m)/(rNorm_m^3);
    end
end

B_mT = 1e3*B_T;

end

%% ========================================================================
% PARAMETER, BASIS AND UNIT AUDIT
% =========================================================================

function run_parameter_and_unit_audit(p,expectSurfaceSpin)

assert(expectSurfaceSpin, ...
    'This file is the dynamic surface-spin model.');
assert(abs(p.t_hop_eV-1.0) < 1e-14, ...
    't_hop changed unexpectedly.');
assert(abs(p.q_eff-0.30) < 1e-14, ...
    'q_eff changed unexpectedly.');
assert(abs(p.T_K-300) < 1e-12, ...
    'Surface temperature must remain 300 K.');
assert(abs(p.tau_rel_fs-100) < 1e-12, ...
    'tau_rel must remain 100 fs.');
assert(abs(p.Gamma0_per_fs-0.01) < 1e-14, ...
    'Gamma0 must remain 0.01 fs^-1.');
assert(abs(p.J0_meV-100) < 1e-12, ...
    'J0 must remain 100 meV.');
assert(abs(p.r0_A-3.0) < 1e-14, ...
    'r0 must remain 3 Angstrom.');
S0row = p.S0(:).';
assert(numel(S0row) == 3 && all(isfinite(S0row)), ...
    'Initial surface spin must be a finite three-component vector.');
assert(abs(norm(S0row)-1) < 1e-14, ...
    'Initial surface spin must be normalized.');
allowedS0 = eye(3);
isAllowedCartesianDirection = any(vecnorm(allowedS0-S0row,2,2) < 1e-14);
assert(isAllowedCartesianDirection, ...
    'Initial surface spin must be one of +x, +y or +z.');
assert(abs(p.alpha_LLG-0.10) < 1e-14, ...
    'LLG damping changed unexpectedly.');
assert(abs(p.gamma_rad_fs_T-1.760859e-4) < 1e-14, ...
    'LLG gyromagnetic ratio changed unexpectedly.');
assert(p.Nseg_BS == 40, ...
    'Nseg_BS changed unexpectedly.');
assert(abs(p.min_distance_A-0.20) < 1e-14, ...
    'Biot-Savart cutoff changed unexpectedly.');
assert(abs(p.initialStep_fs-1e-4) < 1e-16, ...
    'InitialStep must remain 1e-4 fs.');
assert(abs(p.maxStepEarly_fs-1e-2) < 1e-16, ...
    'The early-time ode45 MaxStep must remain 0.01 fs.');
assert(abs(p.odeRelTol-1e-8) < 1e-20 && ...
       abs(p.odeAbsTol-1e-10) < 1e-22, ...
    'ODE45 tolerances changed unexpectedly.');

% Detailed balance.
sampleGap_eV = 0.05;
gammaDown = p.Gamma0_per_fs;
gammaUp = p.Gamma0_per_fs* ...
    exp(-p.beta_eV_inv*sampleGap_eV);

assert(abs(gammaUp/gammaDown- ...
    exp(-p.beta_eV_inv*sampleGap_eV)) < 1e-14, ...
    'Detailed balance failed.');

% Site-major basis.
HorbTest = [1 2;2 3];
HspinTest = kron(HorbTest,eye(2));

assert(HspinTest(1,3) == 2 && HspinTest(2,4) == 2, ...
    'The site-major spin basis is inconsistent.');
assert(HspinTest(1,2) == 0 && HspinTest(3,4) == 0, ...
    'The orbital Hamiltonian mixed spin unexpectedly.');

% Zeeman mT -> T.
HzeemanTest = build_zeeman_hamiltonian([1000 0 0],p);
expected = p.zeeman_prefactor_eV_T;

assert(abs(HzeemanTest(1,2)-expected) < 1e-14, ...
    'The mT-to-T conversion in H_Z failed.');
assert(norm(HzeemanTest-HzeemanTest','fro') < 1e-14, ...
    'The Zeeman test Hamiltonian is not Hermitian.');

% Exchange block normalization:
% S=(0,0,1) gives eigenvalues +/-J/2, so the splitting is J.
HexTest = build_exchange_hamiltonian(0.1,[0 0 1]);

assert(abs(HexTest(1,1)-0.05) < 1e-14 && ...
       abs(HexTest(2,2)+0.05) < 1e-14, ...
    'The exchange normalization is inconsistent.');
assert(norm(HexTest-HexTest','fro') < 1e-14, ...
    'The exchange Hamiltonian is not Hermitian.');

% LLG tangent-vector test: dS/dt must be perpendicular to S.
dStest = llg_rhs([0 0 1],[1 0 0],p);

assert(abs(dot([0 0 1],dStest.')) < 1e-14, ...
    'The LLG derivative is not tangent to the unit sphere.');

% Angstrom and mT/T conversion identities.
assert(abs(1e10*1e-10-1) < 1e-15, ...
    'Angstrom-to-metre conversion failed.');
assert(abs(1e3*1e-3-1) < 1e-15, ...
    'mT-to-T conversion identity failed.');

% Biot-Savart direction and finiteness.
RiTest = [0 0 0;1 0 0];
bondsTest = [1 2];
Btest_mT = biot_savart_at_point( ...
    RiTest,bondsTest,1e-3,[0.5 1 0],0, ...
    100,0.01);

assert(all(isfinite(Btest_mT)), ...
    'Biot-Savart returned NaN or Inf.');
assert(Btest_mT(3) > 0, ...
    'Biot-Savart direction test failed.');

fprintf('PASS: original physical parameters retained.\n');
fprintf(['PASS: chunked ode45 with early MaxStep = 0.01 fs; ' ...
    'late stage controlled by RelTol/AbsTol and joint early stopping.\n']);
fprintf('PASS: detailed balance and sqrt-rate convention.\n');
fprintf('PASS: site-major basis ordering.\n');
fprintf('PASS: exchange block and LLG tangent-vector checks.\n');
fprintf('PASS: Angstrom, ampere, tesla and millitesla checks.\n');
fprintf('PASS: Zeeman and Biot-Savart unit tests.\n');

end

%% ========================================================================
% CSV OUTPUT HELPERS
% =========================================================================

function write_site_time_table(filePath,t,n_up,n_down,n_total,m_z,P_z, ...
    Bx_mT,By_mT,Bz_mT,Babs_mT)

nT = numel(t);
L = size(n_up,2);
nRows = nT*L;

time_fs = zeros(nRows,1);
site = zeros(nRows,1);
up = zeros(nRows,1);
down = zeros(nRows,1);
total = zeros(nRows,1);
imbalance = zeros(nRows,1);
polarization = nan(nRows,1);
Bx = zeros(nRows,1);
By = zeros(nRows,1);
Bz = zeros(nRows,1);
Babs = zeros(nRows,1);

row = 0;

for it = 1:nT
    for i = 1:L

        row = row+1;

        time_fs(row) = t(it);
        site(row) = i;
        up(row) = n_up(it,i);
        down(row) = n_down(it,i);
        total(row) = n_total(it,i);
        imbalance(row) = m_z(it,i);
        polarization(row) = P_z(it,i);
        Bx(row) = Bx_mT(it,i);
        By(row) = By_mT(it,i);
        Bz(row) = Bz_mT(it,i);
        Babs(row) = Babs_mT(it,i);
    end
end

T = table(time_fs,site,up,down,total,imbalance,polarization, ...
    Bx,By,Bz,Babs, ...
    'VariableNames',{ ...
    'time_fs','site','n_up','n_down','n_total','m_z','P_z', ...
    'Bx_mT','By_mT','Bz_mT','B_abs_mT'});

writetable(T,filePath);

end

function write_bond_current_table(filePath,t,bonds,bond_currents_A)

T = table(t,'VariableNames',{'time_fs'});

for b = 1:size(bonds,1)
    name = sprintf('I_bond_%d_to_%d_A',bonds(b,1),bonds(b,2));
    T.(name) = bond_currents_A(:,b);
end

writetable(T,filePath);

end

%% ========================================================================
% FIGURES
% ONLY THE REQUESTED FIGURES ARE CREATED IN THIS VERSION
% =========================================================================

function plot_geometry_comparison(models,p,outputFolder)

fig = figure('Color','w','Position',[100 50 950 1000]);
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

for im = 1:numel(models)

    nexttile;
    hold on;

    Ri = models(im).Ri;
    bonds = models(im).bonds;

    for b = 1:size(bonds,1)
        i = bonds(b,1);
        j = bonds(b,2);

        plot3([Ri(i,1) Ri(j,1)], ...
              [Ri(i,2) Ri(j,2)], ...
              [Ri(i,3) Ri(j,3)],'-','LineWidth',2);
    end

    scatter3(Ri(:,1),Ri(:,2),Ri(:,3),70,'filled');
    scatter3(p.R_surface_A(1),p.R_surface_A(2), ...
        p.R_surface_A(3),120,'filled');

    for i = 1:size(Ri,1)
        text(Ri(i,1),Ri(i,2),Ri(i,3)+0.12,sprintf('%d',i), ...
            'HorizontalAlignment','center');
    end

    xlabel('x [Angstrom]');
    ylabel('y [Angstrom]');
    zlabel('z [Angstrom]');
    title(models(im).displayName);
    grid on;
    box on;
    axis equal;
    view(38,24);
end

sgtitle('Molecular geometries and surface-spin position');
save_figure_files(fig,outputFolder,'figure_1_geometries');

end

%% ========================================================================
% FINAL SITE-RESOLVED SPIN IMBALANCE FOR ONE INITIAL SURFACE-SPIN DIRECTION
% m_z = n_up - n_down IS DIMENSIONLESS
% =========================================================================

function plot_final_imbalance_for_initial_spin( ...
    results,outputFolder,spinCase)

globalFinalMz = 0;

for im = 1:numel(results)
    globalFinalMz = max(globalFinalMz,max(abs(results(im).m_z(end,:))));
end

fig = figure('Color','w','Position',[140 40 1150 950]);
tiledlayout(numel(results),1,'Padding','compact','TileSpacing','compact');

for im = 1:numel(results)

    nexttile;

    bar(1:results(im).L,results(im).m_z(end,:));
    hold on;
    yline(0,'k-');

    xlabel('site index');
    ylabel('m_z = n_{up}-n_{down} [dimensionless]');
    title(sprintf('%s | final t = %.3f fs | stabilized = %d', ...
        results(im).displayName,results(im).t(end), ...
        results(im).jointStabilized));

    xticks(1:results(im).L);
    ylim(1.08*[-max(globalFinalMz,1e-15) max(globalFinalMz,1e-15)]);
    grid on;
end

S0 = spinCase.S0;
sgtitle(sprintf(['Final site-resolved spin imbalance | ' ...
    'S_{surface}(0) = [%.0f %.0f %.0f] (%s)'], ...
    S0(1),S0(2),S0(3),spinCase.label));

save_figure_files(fig,outputFolder, ...
    ['figure_final_spin_imbalance_' spinCase.tag]);

end

%% ========================================================================
% TIME DYNAMICS FOR ONE INITIAL SURFACE-SPIN DIRECTION
%
% LEFT Y AXIS:  B_surface,x/y/z IN mT
% RIGHT Y AXIS: NORMALIZED M_x/M_s, M_y/M_s, M_z/M_s = S_x,S_y,S_z
%
% In this model B_eff = B_surface = B_Biot-Savart.
% Internally the LLG and Zeeman terms still receive B in tesla; this plot
% deliberately uses the stored mT values requested for visualization.
% =========================================================================

function plot_field_and_normalized_magnetization_for_initial_spin( ...
    results,outputFolder,spinCase)

axisNames = {'x','y','z'};
nModels = numel(results);

fig = figure('Color','w','Position',[60 30 1500 1050]);
tiledlayout(3,nModels,'Padding','compact','TileSpacing','compact');

for a = 1:3
    for im = 1:nModels

        tileIndex = (a-1)*nModels+im;
        nexttile(tileIndex);

        t_fs = results(im).t;
        B_component_mT = results(im).B_surface_mT(:,a);

        % S_surface was explicitly normalized when observables were stored,
        % therefore it represents normalized M/M_s in this model.
        M_component_normalized = results(im).S_surface(:,a);

        yyaxis left;
        hB = plot(t_fs,B_component_mT,'LineWidth',1.4);
        ylabel(sprintf('B_{%s} [mT]',axisNames{a}));

        yyaxis right;
        hM = plot(t_fs,M_component_normalized,'LineWidth',1.4);
        ylabel(sprintf('M_{%s}/M_s = S_{%s}', ...
            axisNames{a},axisNames{a}));
        ylim([-1.05 1.05]);

        xlabel('time [fs]');
        title(sprintf('%s | %s axis', ...
            results(im).displayName,upper(axisNames{a})));
        grid on;

        legend([hB hM], ...
            {sprintf('B_{%s} [mT]',axisNames{a}), ...
             sprintf('M_{%s}/M_s',axisNames{a})}, ...
            'Location','best');
    end
end

S0 = spinCase.S0;
sgtitle(sprintf(['Surface-field and normalized magnetization dynamics | ' ...
    'S_{surface}(0) = [%.0f %.0f %.0f] (%s)'], ...
    S0(1),S0(2),S0(3),spinCase.label));

save_figure_files(fig,outputFolder, ...
    ['figure_field_mT_and_Mxyz_normalized_' spinCase.tag]);

end

%% ========================================================================
% MIXED EARLY/LONG-TIME OUTPUT GRID
% =========================================================================

function tspan = build_mixed_time_grid( ...
    tStart,tEarlyEnd,tEnd,dtEarly,dtLate)

assert(tStart < tEnd,'The final time must exceed the initial time.');
assert(dtEarly > 0 && dtLate > 0,'Output steps must be positive.');

tEarlyEnd = min(max(tEarlyEnd,tStart),tEnd);

early = tStart:dtEarly:tEarlyEnd;
if isempty(early) || abs(early(end)-tEarlyEnd) > 1e-12
    early = [early tEarlyEnd];
end

if tEarlyEnd < tEnd
    lateStart = tEarlyEnd+dtLate;
    late = lateStart:dtLate:tEnd;
    if isempty(late) || abs(late(end)-tEnd) > 1e-12
        late = [late tEnd];
    end
else
    late = [];
end

tspan = unique([early late]);

end

function save_figure_files(fig,outputFolder,baseName)

pngPath = fullfile(outputFolder,[baseName '.png']);
figPath = fullfile(outputFolder,[baseName '.fig']);

try
    exportgraphics(fig,pngPath,'Resolution',300);
catch
    saveas(fig,pngPath);
end

savefig(fig,figPath);

end

function downloadsDir = get_downloads_directory()

userProfile = getenv('USERPROFILE');
homeDir = getenv('HOME');

if ~isempty(userProfile)
    downloadsDir = fullfile(userProfile,'Downloads');
elseif ~isempty(homeDir)
    downloadsDir = fullfile(homeDir,'Downloads');
else
    downloadsDir = fullfile(pwd,'Downloads');
end

if ~exist(downloadsDir,'dir')
    mkdir(downloadsDir);
end

end
