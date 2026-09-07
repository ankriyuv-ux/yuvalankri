function soc_zonly_j150()
% HELIX_SOC_BOARD_NN_EXCHANGE_BZONLY_10PS
% =========================================================================
% CLEAN COMPARISON RUN: BOARD-NN SOC + EXCHANGE + SIGNED Z-ONLY LLG
%
% This file is the direct partner of:
%   helix_J0zero_SOC_BOARD_NN_5meV_BZonly_10ps.m
%
% ONLY PHYSICAL CHANGE RELATIVE TO THAT J0=0 CONTROL:
%   J0 = 150 meV and the self-consistent exchange feedback is ON.
%
% EVERYTHING ELSE IS HELD FIXED:
%   - same 13-site helix and contact potential
%   - same initial electron state
%   - same Lindblad bath
%   - same BOARD-NN SOC rule
%       H_SOC,ij = i*(gamma*t)*(S . Rhat_ij),  S=sigma/2
%   - same SOC scan: gamma*t = [0, 0.1, +5, -5] meV
%   - same M(0)=[1 0 0], lambda_LLG=0.6, tau_rel=500 fs
%   - same 10 ps runtime and ode45 controls
%   - no direct electronic Zeeman term
%
% ELECTRON HAMILTONIAN:
%   H_e = H_surface + H_SOC + H_MS[M]
%   H_MS = sum_i J_i M . s_i
%   J_i  = J0 exp(-|R_i-R*|/r0), r0=3 A
%
% EXCHANGE FIELD:
%   B_ex = -(1/(2*mu_surface)) sum_i J_i m_i
%
% IMPORTANT LLG RESTRICTION FOR CLEAN COMPARISON:
%   The electrons evolve with the FULL SOC + exchange Hamiltonian.
%   B_mol and B_ex are both computed in full 3D and saved.
%   But only the SIGNED total z component is passed to the LLG:
%
%       B_LLG = [0, 0, B_eff,z]
%       B_eff,z = B_mol,z + B_ex,z
%
%   Thus x/y are zeroed ONLY at the LLG input, not in H_MS and not in the
%   electronic dynamics. |B_eff,z| is NOT used in the physical dynamics.
%
% PRIMARY QUESTION:
%   At identical BOARD-NN SOC, does turning on J0=150 meV create a
%   significant signed z-directed exchange channel and change Mz relative
%   to the J0=0 control?
%
% NUMERICS:
%   ode45; Refine=1; every accepted solver point retained; no interpolation
%   or downsampling of the trajectory.
% =========================================================================

clc;
close all;
totalTimer = tic;

%% ========================================================================
% 1. PARAMETERS
% =========================================================================

p.t_hop_eV = 1.0;

p.e2_over_4pieps0_eVA = 14.4;
p.q_eff = 0.30;
p.R_contact_A = [0 0 0];

p.Rstar_A = [-1.335757 -0.292878 0.0];

% Localized/surface moment starts exactly along +X.
p.M0 = [1 0 0];
p.M0 = p.M0/norm(p.M0);

p.gamma0_fs_T = 1.760859e-4;
p.lambda = 0.60;

% Exchange held fixed to the strongest previously analysed case so the
% SOC scan can be compared directly with the existing J0=150 meV result.
p.r0_A = 3.0;
p.J0_meV = 150;

% SOC sensitivity scan.  0 meV is the exact baseline control.
p.lambdaSOCList_meV = [0 0.1 5 -5];

p.kB_eV_K = 8.617333262e-5;
p.T_K = 300;
p.beta_eV_inv = 1/(p.kB_eV_K*p.T_K);

p.tau_rel_fs = 500;
p.Gamma0_per_fs = 1/p.tau_rel_fs;

p.hbar_eV_fs = 0.6582119569;
p.muB_eV_T = 5.7883818060e-5;

p.surface_moment_muB = 1.0;
p.surface_moment_eV_T = p.surface_moment_muB*p.muB_eV_T;

% Pilot physical time.
p.tStart_fs = 0;
p.earlyEnd_fs = 50;
p.tEnd_fs = 10000;       % 10 ps
p.chunk_fs = 5000;       % 5 ps chunks after 50 fs

% ode45 controls.
p.RelTol = 1e-8;
p.AbsTol = 1e-10;
p.InitialStep_fs = 1e-4;
p.MaxStepEarly_fs = 1e-2;

% Biot-Savart.
p.Nseg_BS = 40;
p.min_distance_A = 0.20;

p.density_threshold = 1e-8;

% Parallel SOC cases change wall-clock only.
p.useParallelIfAvailable = true;
p.requestedWorkers = min(4,numel(p.lambdaSOCList_meV));

%% ========================================================================
% 2. AUDIT
% =========================================================================

assert(isequal(p.M0,[1 0 0]),'M(0) must be +X.');
assert(abs(p.lambda-0.6) < 1e-14,'lambda must be 0.6.');
assert(abs(p.tau_rel_fs-500) < 1e-12,'tau_rel must be 500 fs.');
assert(any(p.lambdaSOCList_meV==0),'SOC scan must include 0 meV control.');

fprintf('\n============================================================\n');
fprintf('BOARD-NN SOC + EXCHANGE | SIGNED Z-ONLY LLG\n');
fprintf('============================================================\n');
fprintf('J0               = %g meV\n',p.J0_meV);
fprintf('SOC scan          = ');
fprintf('%g ',p.lambdaSOCList_meV);
fprintf('meV\n');
fprintf('M(0)              = [%g %g %g] = +X\n',p.M0);
fprintf('lambda LLG        = %.3f\n',p.lambda);
fprintf('tau_rel           = %.1f fs\n',p.tau_rel_fs);
fprintf('physical time     = %.2f ps\n',p.tEnd_fs/1000);
fprintf('exchange feedback = ON\n');
fprintf('electronic Zeeman = OFF\n');
fprintf('SOC model         = BOARD-NN: i*(gamma*t)*(S . Rhat_ij), S=sigma/2\n');
fprintf('LLG field         = [0,0,B_eff,z] ONLY (SIGNED)\n');
fprintf('solver            = ode45, Refine=1, no downsampling\n');
fprintf('============================================================\n');

%% ========================================================================
% 3. GEOMETRY + BASE STATIC HAMILTONIAN
% =========================================================================

model = build_upright_helix_model();
audit_geometry(model);

% NNN curvature/SOC geometry is built below and audited after construction.

Ri = model.Ri;
bonds = model.bonds;
L = size(Ri,1);
Nspin = 2*L;

model.t_bonds_eV = p.t_hop_eV*model.t_bonds(:);

H_hop = build_hopping_hamiltonian(L,bonds,model.t_bonds_eV);

distance_contact_A = vecnorm(Ri-p.R_contact_A,2,2);
distance_contact_A(distance_contact_A < 1e-12) = 1e-12;

epsilon_site_eV = p.q_eff*p.e2_over_4pieps0_eVA./distance_contact_A;
H_contact = diag(epsilon_site_eV);
H_surface_orbital = H_hop+H_contact;

% Spinful static Hamiltonian without SOC.
model.H_surface_spin_noSOC = kron(H_surface_orbital,eye(2));

% BOARD-NN SOC matrix: unit amplitude in eV is applied later case-by-case.
[model.Hsoc_unit,model.soc_Rhat,model.soc_nn_pairs] = build_board_nn_soc_unit(model);

rhatNorm = vecnorm(model.soc_Rhat,2,2);
fprintf('\nBOARD-NN SOC geometry audit:\n');
fprintf('  NN SOC links        = %d\n',size(model.soc_nn_pairs,1));
fprintf('  |Rhat_ij| range     = %.9f ... %.9f\n',min(rhatNorm),max(rhatNorm));
fprintf('  max |Hsoc_unit|     = %.9f (dimensionless)\n',max(abs(model.Hsoc_unit(:))));

assert(norm(model.Hsoc_unit-model.Hsoc_unit','fro') < 1e-12, ...
    'H_SOC unit matrix is not Hermitian.');

% Precompute indices for local spin observables.
model = precompute_fast_indices(model);

% Biot-Savart kernel at fixed R*.
model.BS_kernel_Rstar = precompute_BS_kernel_at_points( ...
    Ri,bonds,p.Rstar_A,p.Nseg_BS,p.min_distance_A);

%% ========================================================================
% 4. INITIAL ELECTRONIC STATE
% =========================================================================

[rho0,initialDiagnostic] = build_gas_ground_initial_state(H_hop);

fprintf('\nInitial electronic state audit:\n');
fprintf('  trace error          = %.3e\n',initialDiagnostic.traceError);
fprintf('  Hermiticity error    = %.3e\n',initialDiagnostic.hermiticityError);
fprintf('  initial total Mz     = %.3e\n',initialDiagnostic.initialTotalMz);

%% ========================================================================
% 5. SAME FIXED SPIN-INDEPENDENT LINDBLAD BATH FOR ALL SOC CASES
% =========================================================================
% Keeping the bath identical isolates the effect of H_SOC itself.

bath = build_fast_lindblad_bath(H_surface_orbital,p);

% Basic trace-preservation check of the fast bath.
dr0 = apply_fast_lindblad(rho0,bath);
assert(abs(trace(dr0)) < 1e-11,'Fast Lindblad trace check failed.');

%% ========================================================================
% 6. EXCHANGE PROFILE + MATRICES (IDENTICAL IN ALL SOC CASES)
% =========================================================================

J_site_eV = exchange_profile(p.J0_meV,model,p);
exchange = precompute_exchange_matrices(J_site_eV);

%% ========================================================================
% 7. VALIDATE GENERAL BOND CURRENT AT SOC=0
% =========================================================================

case0 = build_soc_case_model(model,0);
Jold = compute_bond_current_scalar_reference(rho0,model);
Jnew = compute_bond_current_general(rho0,case0);

assert(norm(Jold-Jnew) < 1e-13, ...
    'General SOC-capable bond current does not reproduce SOC=0 baseline.');

fprintf('\nPhysical NN backbone-current validation at lambda_SOC=0: PASS\n');

%% ========================================================================
% 8. OUTPUT FOLDER
% =========================================================================

downloadsDir = get_downloads_directory();
runTag = datestr(now,'yyyymmdd_HHMMSS');
outputFolder = fullfile(downloadsDir, ...
    ['SOC_BOARD_NN_EXCHANGE_J0_150_BZONLY_LAMBDA06_10PS_' runTag]);

if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

fprintf('\nOutput folder:\n%s\n',outputFolder);

%% ========================================================================
% 9. PARALLEL AVAILABILITY
% =========================================================================

useParallel = false;

if p.useParallelIfAvailable && license('test','Distrib_Computing_Toolbox')
    try
        pool = gcp('nocreate');
        if isempty(pool)
            pool = parpool('local',p.requestedWorkers);
        end
        useParallel = pool.NumWorkers >= 2;
    catch ME
        fprintf('\nParallel pool unavailable: %s\n',ME.message);
        fprintf('Falling back to serial execution.\n');
        useParallel = false;
    end
end

if useParallel
    fprintf('\nSOC cases will run in PARALLEL.\n');
else
    fprintf('\nSOC cases will run SERIALLY.\n');
end

%% ========================================================================
% 10. RUN SOC CASES
% =========================================================================

nSOC = numel(p.lambdaSOCList_meV);
runs = cell(nSOC,1);

if useParallel
    parfor is = 1:nSOC
        soc_meV = p.lambdaSOCList_meV(is);
        runs{is} = run_one_soc_case( ...
            rho0,soc_meV,J_site_eV,exchange,model,bath,p);
    end
else
    for is = 1:nSOC
        soc_meV = p.lambdaSOCList_meV(is);
        runs{is} = run_one_soc_case( ...
            rho0,soc_meV,J_site_eV,exchange,model,bath,p);
    end
end

%% ========================================================================
% 11. SAVE EACH CASE + BUILD SUMMARY
% =========================================================================

soc_col = zeros(nSOC,1);
Mz_end = zeros(nSOC,1);
maxAbsMz = zeros(nSOC,1);
maxTilt_deg = zeros(nSOC,1);
maxBmolPerp_mT = zeros(nSOC,1);
maxBexPerp_T = zeros(nSOC,1);
maxBeffPerp_T = zeros(nSOC,1);
maxSpinPerp = zeros(nSOC,1);
maxTorque_per_fs = zeros(nSOC,1);
maxAbs_dMz_per_fs = zeros(nSOC,1);
runtime_min = zeros(nSOC,1);
nSolverPoints = zeros(nSOC,1);
BexPerpZ_end_mT = zeros(nSOC,1);
spinPerpZ_end = zeros(nSOC,1);
maxAbsBmolZ_mT = zeros(nSOC,1);
maxAbsBexZ_mT = zeros(nSOC,1);
maxAbsBeffZ_mT = zeros(nSOC,1);
signedBeffZArea_Tfs = zeros(nSOC,1);
absBeffZArea_Tfs = zeros(nSOC,1);
cancelBeffZ_end = zeros(nSOC,1);
Mz_exact_Zonly_end = zeros(nSOC,1);
Mz_abs_headroom_end = zeros(nSOC,1);

for is = 1:nSOC
    E = runs{is};

    soc_col(is) = E.lambdaSOC_meV;
    Mz_end(is) = E.M(end,3);
    maxAbsMz(is) = max(abs(E.M(:,3)));
    maxTilt_deg(is) = max(E.tilt_deg);
    maxBmolPerp_mT(is) = max(E.Bmol_perp_T)*1e3;
    maxBexPerp_T(is) = max(E.Bex_perp_T);
    maxBeffPerp_T(is) = max(E.Beff_perp_T);
    maxSpinPerp(is) = max(E.spin_perp);
    maxTorque_per_fs(is) = max(E.dM_mag_per_fs);
    maxAbs_dMz_per_fs(is) = max(abs(E.dMz_per_fs));
    runtime_min(is) = E.runtime_min;
    nSolverPoints(is) = numel(E.t_fs);
    BexPerpZ_end_mT(is) = 1e3*E.Bex_perp_vec_T(end,3);
    spinPerpZ_end(is) = E.spin_perp_vec(end,3);

    % Z-only observables used by the actual LLG in this run.
    BmolZ_T = 1e-3*E.B_mol_mT(:,3);
    BexZ_T = E.B_ex_T(:,3);
    BeffZ_T = E.B_eff_T(:,3);
    Az_signed = cumtrapz(E.t_fs(:),BeffZ_T);
    Az_abs = cumtrapz(E.t_fs(:),abs(BeffZ_T));

    E.cumBeffZ_signed_Tfs = Az_signed;
    E.cumBeffZ_abs_Tfs = Az_abs;
    E.BeffZ_cancellationRatio = abs(Az_signed)./max(Az_abs,eps);
    E.BeffZ_cancellationRatio(1) = NaN;
    E.Mz_exact_Zonly = tanh(p.lambda*p.gamma0_fs_T*Az_signed);
    E.Mz_abs_headroom = tanh(p.lambda*p.gamma0_fs_T*Az_abs);
    runs{is} = E;

    maxAbsBmolZ_mT(is) = 1e3*max(abs(BmolZ_T));
    maxAbsBexZ_mT(is) = 1e3*max(abs(BexZ_T));
    maxAbsBeffZ_mT(is) = 1e3*max(abs(BeffZ_T));
    signedBeffZArea_Tfs(is) = Az_signed(end);
    absBeffZArea_Tfs(is) = Az_abs(end);
    cancelBeffZ_end(is) = abs(Az_signed(end))/max(Az_abs(end),eps);
    Mz_exact_Zonly_end(is) = E.Mz_exact_Zonly(end);
    Mz_abs_headroom_end(is) = E.Mz_abs_headroom(end);

    tagS = safe_number_tag(E.lambdaSOC_meV);

    save(fullfile(outputFolder, ...
        sprintf('SOC_%s_meV_COMPACT.mat',tagS)),'E','-v7.3');

    write_case_csv(E,outputFolder);
end

i0 = find(abs(soc_col) < 1e-15,1,'first');
if isempty(i0)
    error('SOC=0 control missing from completed runs.');
end
Mz_delta_from_zero = Mz_end - Mz_end(i0);
SOC_to_hopping = (1e-3*soc_col)/p.t_hop_eV;

summaryTable = table( ...
    soc_col,SOC_to_hopping,Mz_end,Mz_delta_from_zero,maxAbsMz,maxTilt_deg, ...
    maxBmolPerp_mT,maxBexPerp_T,maxBeffPerp_T,maxSpinPerp, ...
    BexPerpZ_end_mT,spinPerpZ_end, ...
    maxAbsBmolZ_mT,maxAbsBexZ_mT,maxAbsBeffZ_mT, ...
    signedBeffZArea_Tfs,absBeffZArea_Tfs,cancelBeffZ_end, ...
    Mz_exact_Zonly_end,Mz_abs_headroom_end, ...
    maxTorque_per_fs,maxAbs_dMz_per_fs,nSolverPoints,runtime_min, ...
    'VariableNames',{ ...
    'lambda_SOC_meV','lambdaSOC_over_t','Mz_end','delta_Mz_end_vs_SOC0', ...
    'max_abs_Mz','max_tilt_deg', ...
    'max_Bmol_perp_mT','max_Bex_perp_T','max_Beff_perp_T', ...
    'max_electronic_spin_perp','Bex_perp_z_end_mT','spin_perp_z_end', ...
    'max_abs_Bmol_z_mT','max_abs_Bex_z_mT','max_abs_Beff_z_mT', ...
    'signed_Beff_z_area_Tfs','abs_Beff_z_area_Tfs','Beff_z_cancellation_ratio_end', ...
    'Mz_exact_Zonly_end','Mz_abs_headroom_end', ...
    'max_dM_dt_per_fs','max_abs_dMz_dt_per_fs', ...
    'n_solver_points','runtime_min'});

writetable(summaryTable,fullfile(outputFolder,'BOARD_NN_EXCHANGE_BZONLY_SUMMARY.csv'));

save(fullfile(outputFolder,'COMPLETE_BOARD_NN_EXCHANGE_BZONLY.mat'), ...
    'runs','summaryTable','p','model','J_site_eV','initialDiagnostic','-v7.3');

%% ========================================================================
% 12. PRINT NUMERICAL SUMMARY
% =========================================================================

fprintf('\n============================================================\n');
fprintf('BOARD-NN SOC + EXCHANGE | Z-ONLY SUMMARY\n');
fprintf('============================================================\n');
disp(summaryTable);

% Enhancement relative to SOC=0 baseline.
if ~isempty(i0)
    baseBex = maxBexPerp_T(i0);
    baseBeff = maxBeffPerp_T(i0);
    baseMz = maxAbsMz(i0);

    fprintf('\nEnhancement relative to SOC=0:\n');
    for is = 1:nSOC
        fprintf(['  SOC=%g meV: Bex_perp x %.4g | Beff_perp x %.4g | ' ...
                 'max|Mz| x %.4g\n'], ...
            soc_col(is), ...
            maxBexPerp_T(is)/max(baseBex,eps), ...
            maxBeffPerp_T(is)/max(baseBeff,eps), ...
            maxAbsMz(is)/max(baseMz,eps));
    end
end

fprintf('\nSOC sign / effective chirality control (final-time changes vs SOC=0):\n');
absLevels = unique(abs(soc_col(soc_col~=0)));
for ia = 1:numel(absLevels)
    a = absLevels(ia);
    ip = find(abs(soc_col-a) < 1e-12,1,'first');
    im = find(abs(soc_col+a) < 1e-12,1,'first');
    if ~isempty(ip) && ~isempty(im)
        dp = Mz_delta_from_zero(ip);
        dm = Mz_delta_from_zero(im);
        fprintf('  |SOC|=%g meV: deltaMz(+)= %+.6e | deltaMz(-)= %+.6e',a,dp,dm);
        if dp*dm < 0
            fprintf('  -> OPPOSITE SIGNS (chirality-odd response present)\n');
        else
            fprintf('  -> not opposite at final time\n');
        end
    end
end

%% ========================================================================
% 13. FIGURES
% =========================================================================

labels = arrayfun(@(x)sprintf('SOC = %g meV',x),soc_col, ...
    'UniformOutput',false);

% ---- FIGURE 0: helix + BOARD-NN SOC bond directions ---------------------
fig0 = figure('Color','w','Position',[60 60 1200 850]);
plot3(model.Ri(:,1),model.Ri(:,2),model.Ri(:,3),'-o','LineWidth',1.2);
hold on;
for ii = 1:size(model.soc_nn_pairs,1)
    i = model.soc_nn_pairs(ii,1); j = model.soc_nn_pairs(ii,2);
    rmid = 0.5*(model.Ri(i,:)+model.Ri(j,:));
    rh = model.soc_Rhat(ii,:);
    quiver3(rmid(1),rmid(2),rmid(3),rh(1),rh(2),rh(3),0.55,'LineWidth',1.0);
end
grid on; axis equal;
xlabel('x [A]'); ylabel('y [A]'); zlabel('z [A]');
title('Helix geometry and BOARD-NN SOC bond directions Rhat_{ij}');
view(35,25);
save_figure(fig0,outputFolder,'FIG00_BOARD_NN_SOC_GEOMETRY.png');

% ---- FIGURE 1: Mz(t) ----------------------------------------------------
fig1 = figure('Color','w','Position',[80 80 1400 760]);
hold on;
for is = 1:nSOC
    plot(runs{is}.t_fs/1000,runs{is}.M(:,3),'LineWidth',1.35, ...
        'DisplayName',labels{is});
end
grid on;
xlabel('time [ps]');
ylabel('M_z');
title(sprintf('BOARD-NN SOC: M_z(t) | J_0=%g meV | signed B_{eff,z}-only LLG',p.J0_meV));
legend('Location','best');
xlim([0 p.tEnd_fs/1000]);
save_figure(fig1,outputFolder,'FIG01_Mz_SOC_SCAN.png');

% ---- FIGURE 2: Bex perpendicular ----------------------------------------
fig2 = figure('Color','w','Position',[100 100 1400 760]);
hold on;
for is = 1:nSOC
    plot(runs{is}.t_fs/1000,1e3*runs{is}.Bex_perp_T,'LineWidth',1.35, ...
        'DisplayName',labels{is});
end
grid on;
xlabel('time [ps]');
ylabel('|B_{ex,\perp}| [mT]');
title('Transverse exchange field diagnostic (full vector saved)');
legend('Location','best');
xlim([0 p.tEnd_fs/1000]);
save_figure(fig2,outputFolder,'FIG02_Bex_PERP_SOC_SCAN.png');

% ---- FIGURE 3: Beff perpendicular ---------------------------------------
fig3 = figure('Color','w','Position',[120 120 1400 760]);
hold on;
for is = 1:nSOC
    plot(runs{is}.t_fs/1000,1e3*runs{is}.Beff_perp_T,'LineWidth',1.35, ...
        'DisplayName',labels{is});
end
grid on;
xlabel('time [ps]');
ylabel('|B_{eff,\perp}| [mT]');
title('Full B_{eff,\perp} diagnostic (LLG uses z only)');
legend('Location','best');
xlim([0 p.tEnd_fs/1000]);
save_figure(fig3,outputFolder,'FIG03_Beff_PERP_SOC_SCAN.png');

% ---- FIGURE 4: electronic transverse spin -------------------------------
fig4 = figure('Color','w','Position',[140 140 1400 760]);
hold on;
for is = 1:nSOC
    plot(runs{is}.t_fs/1000,runs{is}.spin_perp,'LineWidth',1.35, ...
        'DisplayName',labels{is});
end
grid on;
xlabel('time [ps]');
ylabel('|m_{e,\perp}|');
title('Electronic spin component transverse to M');
legend('Location','best');
xlim([0 p.tEnd_fs/1000]);
save_figure(fig4,outputFolder,'FIG04_ELECTRON_SPIN_PERP_SOC_SCAN.png');

% ---- FIGURE 5: LLG torque magnitude ------------------------------------
fig5 = figure('Color','w','Position',[160 160 1400 760]);
hold on;
for is = 1:nSOC
    plot(runs{is}.t_fs/1000,runs{is}.dM_mag_per_fs,'LineWidth',1.35, ...
        'DisplayName',labels{is});
end
grid on;
xlabel('time [ps]');
ylabel('|dM/dt| [fs^{-1}]');
title('Actual LLG rotation rate / torque diagnostic');
legend('Location','best');
xlim([0 p.tEnd_fs/1000]);
save_figure(fig5,outputFolder,'FIG05_LLG_TORQUE_SOC_SCAN.png');

% ---- FIGURE 6: summary vs SOC ------------------------------------------
[socPlot,iSocPlot] = sort(soc_col);
fig6 = figure('Color','w','Position',[180 100 1450 900]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile;
plot(socPlot,1e3*maxBexPerp_T(iSocPlot),'-o','LineWidth',1.4);
grid on;
xlabel('\lambda_{SOC} [meV]');
ylabel('max |B_{ex,\perp}| [mT]');
title('Transverse exchange field');

nexttile;
plot(socPlot,1e3*maxBeffPerp_T(iSocPlot),'-o','LineWidth',1.4);
grid on;
xlabel('\lambda_{SOC} [meV]');
ylabel('max |B_{eff,\perp}| [mT]');
title('Full B_{eff,\perp} diagnostic (LLG uses z only)');

nexttile;
plot(socPlot,maxAbsMz(iSocPlot),'-o','LineWidth',1.4);
grid on;
xlabel('\lambda_{SOC} [meV]');
ylabel('max |M_z|');
title('Out-of-plane moment response');

nexttile;
plot(socPlot,maxTilt_deg(iSocPlot),'-o','LineWidth',1.4);
grid on;
xlabel('\lambda_{SOC} [meV]');
ylabel('max tilt from +X [deg]');
title('Surface-moment reorientation');

sgtitle(sprintf('BOARD-NN SOC summary | J_0=%g meV | Z-only LLG | \lambda_{LLG}=%.1f', ...
    p.J0_meV,p.lambda));
save_figure(fig6,outputFolder,'FIG06_SOC_SUMMARY.png');

% ---- FIGURE 7: final chirality/sign-control summary ----------------------
fig7 = figure('Color','w','Position',[200 100 1400 900]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

nexttile;
plot(socPlot,Mz_delta_from_zero(iSocPlot),'-o','LineWidth',1.4);
hold on; yline(0,'--');
grid on;
xlabel('\lambda_{SOC} [meV]');
ylabel('\Delta M_z(end)');
title('Final M_z change relative to SOC=0');

nexttile;
plot(socPlot,BexPerpZ_end_mT(iSocPlot),'-o','LineWidth',1.4);
hold on; yline(0,'--');
grid on;
xlabel('\lambda_{SOC} [meV]');
ylabel('B_{ex,\perp,z}(end) [mT]');
title('Signed transverse exchange-field z component');

nexttile;
plot(socPlot,spinPerpZ_end(iSocPlot),'-o','LineWidth',1.4);
hold on; yline(0,'--');
grid on;
xlabel('\lambda_{SOC} [meV]');
ylabel('s_{e,\perp,z}(end)');
title('Signed transverse electronic-spin z component');

sgtitle('SOC sign control: +\lambda versus -\lambda');
save_figure(fig7,outputFolder,'FIG07_SOC_SIGN_CHIRALITY_CONTROL.png');

% ---- FIGURE 8: all M components for BOARD value +/-5 meV ----------------
iPlus = find(abs(soc_col-5)<1e-12,1,'first');
iMinus = find(abs(soc_col+5)<1e-12,1,'first');
Ep = runs{iPlus};
Em = runs{iMinus};
fig8 = figure('Color','w','Position',[220 80 1500 900]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

nexttile; plot(Ep.t_fs/1000,Ep.M(:,1),'LineWidth',1.3); grid on; ylabel('M_x'); title('+5 meV'); xlim([0 p.tEnd_fs/1000]);
nexttile; plot(Ep.t_fs/1000,Ep.M(:,2),'LineWidth',1.3); grid on; ylabel('M_y'); title('+5 meV'); xlim([0 p.tEnd_fs/1000]);
nexttile; plot(Ep.t_fs/1000,Ep.M(:,3),'LineWidth',1.3); grid on; ylabel('M_z'); title('+5 meV'); xlim([0 p.tEnd_fs/1000]);

nexttile; plot(Em.t_fs/1000,Em.M(:,1),'LineWidth',1.3); grid on; ylabel('M_x'); xlabel('time [ps]'); title('-5 meV'); xlim([0 p.tEnd_fs/1000]);
nexttile; plot(Em.t_fs/1000,Em.M(:,2),'LineWidth',1.3); grid on; ylabel('M_y'); xlabel('time [ps]'); title('-5 meV'); xlim([0 p.tEnd_fs/1000]);
nexttile; plot(Em.t_fs/1000,Em.M(:,3),'LineWidth',1.3); grid on; ylabel('M_z'); xlabel('time [ps]'); title('-5 meV'); xlim([0 p.tEnd_fs/1000]);

sgtitle('Moment components for central SOC sign pair | M(0)=+X');
save_figure(fig8,outputFolder,'FIG08_M_COMPONENTS_PLUS_MINUS_5meV.png');

% ---- FIGURE 9: field decomposition for BOARD value +5 meV ----------------
Es = Ep;
fig9 = figure('Color','w','Position',[240 100 1450 920]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

nexttile;
plot(Es.t_fs/1000,Es.B_mol_mT(:,1),'LineWidth',1.1,'DisplayName','B_{mol,x}');
hold on;
plot(Es.t_fs/1000,Es.B_mol_mT(:,2),'LineWidth',1.1,'DisplayName','B_{mol,y}');
plot(Es.t_fs/1000,Es.B_mol_mT(:,3),'LineWidth',1.1,'DisplayName','B_{mol,z}');
grid on; ylabel('B_{mol} [mT]'); legend('Location','best');
xlim([0 p.tEnd_fs/1000]);
title('Molecular Biot-Savart field');

nexttile;
plot(Es.t_fs/1000,Es.B_ex_T(:,1),'LineWidth',1.1,'DisplayName','B_{ex,x}');
hold on;
plot(Es.t_fs/1000,Es.B_ex_T(:,2),'LineWidth',1.1,'DisplayName','B_{ex,y}');
plot(Es.t_fs/1000,Es.B_ex_T(:,3),'LineWidth',1.1,'DisplayName','B_{ex,z}');
grid on; ylabel('B_{ex} [T]'); legend('Location','best');
xlim([0 p.tEnd_fs/1000]);
title('Exchange effective field');

nexttile;
plot(Es.t_fs/1000,1e3*Es.Bex_perp_T,'LineWidth',1.2, ...
    'DisplayName','|B_{ex,\perp}|');
hold on;
plot(Es.t_fs/1000,1e3*Es.Beff_perp_T,'LineWidth',1.2, ...
    'DisplayName','|B_{eff,\perp}|');
grid on; ylabel('transverse field [mT]'); xlabel('time [ps]');
legend('Location','best'); xlim([0 p.tEnd_fs/1000]);
title('Full transverse fields (diagnostic; LLG uses only B_{eff,z})');

sgtitle('Field decomposition | BOARD-NN SOC = +5 meV | Z-only LLG');
save_figure(fig9,outputFolder,'FIG09_FIELD_DECOMPOSITION_PLUS_5meV.png');

% ---- FIGURE 10: the actual Z field entering LLG -------------------------
fig10 = figure('Color','w','Position',[260 70 1500 930]);
tiledlayout(nSOC,1,'TileSpacing','compact','Padding','compact');
for is=1:nSOC
    E=runs{is}; nexttile;
    plot(E.t_fs/1000,E.B_mol_mT(:,3),'LineWidth',1.15,'DisplayName','B_{mol,z}'); hold on;
    plot(E.t_fs/1000,1e3*E.B_ex_T(:,3),'LineWidth',1.15,'DisplayName','B_{ex,z}');
    plot(E.t_fs/1000,1e3*E.B_eff_T(:,3),'LineWidth',1.25,'DisplayName','B_{eff,z}');
    yline(0,'--'); grid on; ylabel('B_z [mT]'); title(labels{is});
    legend('Location','best');
    if is==nSOC, xlabel('time [ps]'); end
end
sgtitle('Actual signed Z channel: B_{eff,z}=B_{mol,z}+B_{ex,z}');
save_figure(fig10,outputFolder,'FIG10_Z_FIELD_DECOMPOSITION_ALL_SOC.png');

% ---- FIGURE 11: accumulated signed vs absolute B_eff,z ------------------
fig11 = figure('Color','w','Position',[280 70 1500 930]);
tiledlayout(nSOC,1,'TileSpacing','compact','Padding','compact');
for is=1:nSOC
    E=runs{is}; nexttile;
    plot(E.t_fs/1000,E.cumBeffZ_signed_Tfs,'LineWidth',1.25,'DisplayName','int B_{eff,z} dt'); hold on;
    plot(E.t_fs/1000,E.cumBeffZ_abs_Tfs,'LineWidth',1.25,'DisplayName','int |B_{eff,z}| dt');
    grid on; ylabel('T fs'); title(labels{is}); legend('Location','best');
    if is==nSOC, xlabel('time [ps]'); end
end
sgtitle('Signed versus absolute accumulated total Z field');
save_figure(fig11,outputFolder,'FIG11_ACCUMULATED_BEFF_Z.png');

% ---- FIGURE 12: cancellation ratio for actual B_eff,z -------------------
fig12 = figure('Color','w','Position',[300 100 1400 760]); hold on;
for is=1:nSOC
    E=runs{is};
    plot(E.t_fs/1000,E.BeffZ_cancellationRatio,'LineWidth',1.25,'DisplayName',labels{is});
end
grid on; xlabel('time [ps]'); ylabel('|int B_{eff,z}dt| / int|B_{eff,z}|dt');
title('Temporal cancellation of the actual total Z field entering LLG');
legend('Location','best'); xlim([0 p.tEnd_fs/1000]);
save_figure(fig12,outputFolder,'FIG12_BEFF_Z_CANCELLATION_RATIO.png');

% ---- FIGURE 13: exact Z-only consistency + |Bz| headroom ----------------
fig13 = figure('Color','w','Position',[320 80 1500 900]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile; hold on;
for is=1:nSOC
    E=runs{is};
    plot(E.t_fs/1000,E.M(:,3),'LineWidth',1.25,'DisplayName',[labels{is} ' LLG']);
    plot(E.t_fs/1000,E.Mz_exact_Zonly,'--','LineWidth',1.0,'DisplayName',[labels{is} ' integral']);
end
grid on; ylabel('M_z'); title('Actual Z-only LLG versus accumulated signed B_{eff,z}'); legend('Location','best');
nexttile; hold on;
for is=1:nSOC
    E=runs{is};
    plot(E.t_fs/1000,E.M(:,3),'LineWidth',1.15,'DisplayName',[labels{is} ' signed']);
    plot(E.t_fs/1000,E.Mz_abs_headroom,'--','LineWidth',1.0,'DisplayName',[labels{is} ' |B_{eff,z}| headroom']);
end
grid on; xlabel('time [ps]'); ylabel('M_z');
title('|B_{eff,z}| is diagnostic headroom ONLY; not used in dynamics'); legend('Location','best');
save_figure(fig13,outputFolder,'FIG13_ZONLY_CONSISTENCY_AND_HEADROOM.png');

%% ========================================================================
% 14. MODEL NOTE
% =========================================================================

fid = fopen(fullfile(outputFolder,'MODEL_NOTE_SOC.txt'),'w');
if fid >= 0
    fprintf(fid,'BOARD-NN SOC + EXCHANGE + SIGNED Z-ONLY LLG\n\n');
    fprintf(fid,'Direct comparison partner of the previous J0=0 BOARD-NN run.\n');
    fprintf(fid,'Only J0 is changed to 150 meV; all other numerical/model settings are held fixed.\n');
    fprintf(fid,'H_SOC,ij = i*(gamma*t)*(S.Rhat_ij), S=sigma/2, gamma*t = 0,0.1,+5,-5 meV.\n');
    fprintf(fid,'H_MS[M] is fully active in the electronic Hamiltonian.\n');
    fprintf(fid,'B_mol and B_ex are computed in full 3D and saved.\n');
    fprintf(fid,'LLG receives only signed [0,0,B_eff,z], where B_eff=B_mol+B_ex.\n');
    fprintf(fid,'No direct electronic Zeeman term. |B_eff,z| is diagnostic only.\n');
    fclose(fid);
end

fprintf('\n============================================================\n');
fprintf('BOARD-NN SOC + EXCHANGE Z-ONLY RUN COMPLETE\n');
fprintf('Output folder:\n%s\n',outputFolder);
fprintf('Total wall time = %.2f min\n',toc(totalTimer)/60);
fprintf('============================================================\n');

end


%% =========================================================================
% ONE SOC CASE
% =========================================================================

function E = run_one_soc_case(rho0,soc_meV,J_site_eV,exchange,model,bath,p)

caseTimer = tic;

fprintf('\n------------------------------------------------------------\n');
fprintf('START SOC = %g meV\n',soc_meV);
fprintf('------------------------------------------------------------\n');

caseModel = build_soc_case_model(model,soc_meV);

Nspin = size(caseModel.H_static_spin,1);
Nrho = Nspin^2;

rhs = @(~,y) coupled_soc_exchange_rhs( ...
    y,J_site_eV,exchange,caseModel,bath,p);

y0 = [rho0(:);p.M0(:)];

% Accumulated observables -- EVERY accepted ode45 point.
tObs = [];
MObs = [];
BmolObs = [];
BexObs = [];
BeffObs = [];
spinObs = [];
BmolPerpObs = [];
BexPerpObs = [];
BeffPerpObs = [];
spinPerpObs = [];
dMmagObs = [];
dMzObs = [];
tiltObs = [];
traceObs = [];
hermObs = [];
MnormObs = [];

% PHASE A: 0-50 fs conservative.
optsEarly = odeset( ...
    'InitialStep',p.InitialStep_fs, ...
    'MaxStep',p.MaxStepEarly_fs, ...
    'RelTol',p.RelTol, ...
    'AbsTol',p.AbsTol, ...
    'Refine',1, ...
    'Stats','off');

[tSeg,YSeg] = ode45(rhs,[p.tStart_fs p.earlyEnd_fs],y0,optsEarly);

[tObs,MObs,BmolObs,BexObs,BeffObs,spinObs, ...
 BmolPerpObs,BexPerpObs,BeffPerpObs,spinPerpObs, ...
 dMmagObs,dMzObs,tiltObs,traceObs,hermObs,MnormObs,lastState] = ...
    append_soc_segment(tSeg,YSeg,false, ...
    tObs,MObs,BmolObs,BexObs,BeffObs,spinObs, ...
    BmolPerpObs,BexPerpObs,BeffPerpObs,spinPerpObs, ...
    dMmagObs,dMzObs,tiltObs,traceObs,hermObs,MnormObs, ...
    J_site_eV,caseModel,p);

fprintf('  SOC=%g meV: reached %.3f ps | points=%d\n', ...
    soc_meV,p.earlyEnd_fs/1000,numel(tSeg));

% PHASE B: >50 fs adaptive, no user MaxStep.
optsLate = odeset( ...
    'RelTol',p.RelTol, ...
    'AbsTol',p.AbsTol, ...
    'Refine',1, ...
    'Stats','off');

currentTime = p.earlyEnd_fs;

while currentTime < p.tEnd_fs-1e-12
    nextEnd = min(currentTime+p.chunk_fs,p.tEnd_fs);
    segTimer = tic;

    [tSeg,YSeg] = ode45(rhs,[currentTime nextEnd],lastState,optsLate);

    [tObs,MObs,BmolObs,BexObs,BeffObs,spinObs, ...
     BmolPerpObs,BexPerpObs,BeffPerpObs,spinPerpObs, ...
     dMmagObs,dMzObs,tiltObs,traceObs,hermObs,MnormObs,lastState] = ...
        append_soc_segment(tSeg,YSeg,true, ...
        tObs,MObs,BmolObs,BexObs,BeffObs,spinObs, ...
        BmolPerpObs,BexPerpObs,BeffPerpObs,spinPerpObs, ...
        dMmagObs,dMzObs,tiltObs,traceObs,hermObs,MnormObs, ...
        J_site_eV,caseModel,p);

    currentTime = nextEnd;

    dt = diff(tSeg);
    if isempty(dt)
        dtMed = NaN;
    else
        dtMed = median(dt);
    end

    fprintf('  SOC=%g meV: reached %5.2f ps | points=%7d | median dt=%g fs | %.2f min\n', ...
        soc_meV,currentTime/1000,numel(tSeg),dtMed,toc(segTimer)/60);
end

rhoFinal = reshape(lastState(1:Nrho),Nspin,Nspin);
rhoFinal = 0.5*(rhoFinal+rhoFinal');

E.lambdaSOC_meV = soc_meV;
E.t_fs = tObs;
E.M = MObs;
E.B_mol_mT = BmolObs;
E.B_ex_T = BexObs;
E.B_eff_T = BeffObs;
E.B_LLG_T = [zeros(size(BeffObs,1),2), BeffObs(:,3)];
E.totalElectronicSpin = spinObs;
E.Bmol_perp_T = BmolPerpObs;
E.Bex_perp_T = BexPerpObs;
E.Beff_perp_T = BeffPerpObs;
E.spin_perp = spinPerpObs;

% Signed transverse vectors (not only magnitudes), useful for +/-SOC test.
E.Bmol_perp_vec_T = (1e-3*BmolObs) - sum((1e-3*BmolObs).*MObs,2).*MObs;
E.Bex_perp_vec_T  = BexObs - sum(BexObs.*MObs,2).*MObs;
E.Beff_perp_vec_T = BeffObs - sum(BeffObs.*MObs,2).*MObs;
E.spin_perp_vec    = spinObs - sum(spinObs.*MObs,2).*MObs;

E.dM_mag_per_fs = dMmagObs;
E.dMz_per_fs = dMzObs;
E.tilt_deg = tiltObs;
E.traceError = traceObs;
E.hermiticityError = hermObs;
E.MnormError = MnormObs;
E.finalRho = rhoFinal;
E.finalObs = density_observables_from_rho(rhoFinal,p);
E.runtime_min = toc(caseTimer)/60;

fprintf('COMPLETE SOC=%g meV | runtime %.2f min | final Mz=%+.6e | max Bex_perp=%.6e T\n', ...
    soc_meV,E.runtime_min,E.M(end,3),max(E.Bex_perp_T));

end


%% =========================================================================
% COUPLED RHS
% =========================================================================

function dy = coupled_soc_exchange_rhs(y,J_site_eV,exchange,model,bath,p)

Nspin = size(model.H_static_spin,1);
Nrho = Nspin^2;

rho = reshape(y(1:Nrho),Nspin,Nspin);
rho = 0.5*(rho+rho');

M = real(y(Nrho+1:Nrho+3));
M = M(:).';
M = M/max(norm(M),1e-15);

% General charge current including spin-dependent SOC hopping.
J_eV = compute_bond_current_general(rho,model);

B_mol_mT = apply_BS_kernel(model.BS_kernel_Rstar,J_eV);
B_mol_T = 1e-3*B_mol_mT;

H_exchange = M(1)*exchange.Hx + M(2)*exchange.Hy + M(3)*exchange.Hz;

% NO direct electronic Zeeman.
H_total = model.H_static_spin + H_exchange;

drho_H = -1i/p.hbar_eV_fs*(H_total*rho-rho*H_total);
drho_D = apply_fast_lindblad(rho,bath);
drho = drho_H+drho_D;

B_ex_T = compute_Bex_fast(rho,J_site_eV,model,p);
B_eff_T = B_mol_T+B_ex_T;

% CLEAN COMPARISON RESTRICTION: only SIGNED total z component drives LLG.
% x/y remain fully present in H_MS and in the saved electronic/exchange data.
B_LLG_T = [0 0 B_eff_T(3)];
precession = cross(M,B_LLG_T);
damping = cross(M,cross(M,B_LLG_T));
dM = p.gamma0_fs_T*(precession-p.lambda*damping);

dy = [drho(:);dM(:)];

end


%% =========================================================================
% APPEND ALL ACCEPTED SOLVER POINTS
% =========================================================================

function [tObs,MObs,BmolObs,BexObs,BeffObs,spinObs, ...
          BmolPerpObs,BexPerpObs,BeffPerpObs,spinPerpObs, ...
          dMmagObs,dMzObs,tiltObs,traceObs,hermObs,MnormObs,lastState] = ...
    append_soc_segment(tSeg,YSeg,dropFirst, ...
    tObs,MObs,BmolObs,BexObs,BeffObs,spinObs, ...
    BmolPerpObs,BexPerpObs,BeffPerpObs,spinPerpObs, ...
    dMmagObs,dMzObs,tiltObs,traceObs,hermObs,MnormObs, ...
    J_site_eV,model,p)

if dropFirst
    tSeg = tSeg(2:end);
    YSeg = YSeg(2:end,:);
end

n = numel(tSeg);
Nspin = size(model.H_static_spin,1);
Nrho = Nspin^2;

Mseg = zeros(n,3);
BmolSeg = zeros(n,3);
BexSeg = zeros(n,3);
BeffSeg = zeros(n,3);
spinSeg = zeros(n,3);
BmolPerpSeg = zeros(n,1);
BexPerpSeg = zeros(n,1);
BeffPerpSeg = zeros(n,1);
spinPerpSeg = zeros(n,1);
dMmagSeg = zeros(n,1);
dMzSeg = zeros(n,1);
tiltSeg = zeros(n,1);
traceSeg = zeros(n,1);
hermSeg = zeros(n,1);
MnormSeg = zeros(n,1);

for it = 1:n
    rhoRaw = reshape(YSeg(it,1:Nrho).',Nspin,Nspin);
    rho = 0.5*(rhoRaw+rhoRaw');

    Mraw = real(YSeg(it,Nrho+1:Nrho+3));
    Mraw = Mraw(:).';
    Mnorm = norm(Mraw);
    M = Mraw/max(Mnorm,1e-15);

    J_eV = compute_bond_current_general(rho,model);
    Bmol_mT = apply_BS_kernel(model.BS_kernel_Rstar,J_eV);
    Bmol_T = 1e-3*Bmol_mT;
    Bex_T = compute_Bex_fast(rho,J_site_eV,model,p);
    Beff_T = Bmol_T+Bex_T;

    obs = density_observables_from_rho(rho,p);
    sTot = [sum(obs.m_x),sum(obs.m_y),sum(obs.m_z)];

    BmolPerpVec = Bmol_T-dot(Bmol_T,M)*M;
    BexPerpVec = Bex_T-dot(Bex_T,M)*M;
    BeffPerpVec = Beff_T-dot(Beff_T,M)*M;
    spinPerpVec = sTot-dot(sTot,M)*M;

    B_LLG_T = [0 0 Beff_T(3)];
    precession = cross(M,B_LLG_T);
    damping = cross(M,cross(M,B_LLG_T));
    dM = p.gamma0_fs_T*(precession-p.lambda*damping);

    Mseg(it,:) = M;
    BmolSeg(it,:) = Bmol_mT;
    BexSeg(it,:) = Bex_T;
    BeffSeg(it,:) = Beff_T;
    spinSeg(it,:) = sTot;
    BmolPerpSeg(it) = norm(BmolPerpVec);
    BexPerpSeg(it) = norm(BexPerpVec);
    BeffPerpSeg(it) = norm(BeffPerpVec);
    spinPerpSeg(it) = norm(spinPerpVec);
    dMmagSeg(it) = norm(dM);
    dMzSeg(it) = dM(3);
    tiltSeg(it) = acosd(max(-1,min(1,M(1))));
    traceSeg(it) = abs(real(trace(rho))-1);
    hermSeg(it) = norm(rhoRaw-rhoRaw','fro');
    MnormSeg(it) = abs(Mnorm-1);
end

tObs = [tObs;tSeg]; %#ok<AGROW>
MObs = [MObs;Mseg]; %#ok<AGROW>
BmolObs = [BmolObs;BmolSeg]; %#ok<AGROW>
BexObs = [BexObs;BexSeg]; %#ok<AGROW>
BeffObs = [BeffObs;BeffSeg]; %#ok<AGROW>
spinObs = [spinObs;spinSeg]; %#ok<AGROW>
BmolPerpObs = [BmolPerpObs;BmolPerpSeg]; %#ok<AGROW>
BexPerpObs = [BexPerpObs;BexPerpSeg]; %#ok<AGROW>
BeffPerpObs = [BeffPerpObs;BeffPerpSeg]; %#ok<AGROW>
spinPerpObs = [spinPerpObs;spinPerpSeg]; %#ok<AGROW>
dMmagObs = [dMmagObs;dMmagSeg]; %#ok<AGROW>
dMzObs = [dMzObs;dMzSeg]; %#ok<AGROW>
tiltObs = [tiltObs;tiltSeg]; %#ok<AGROW>
traceObs = [traceObs;traceSeg]; %#ok<AGROW>
hermObs = [hermObs;hermSeg]; %#ok<AGROW>
MnormObs = [MnormObs;MnormSeg]; %#ok<AGROW>

lastState = YSeg(end,:).';

end


%% =========================================================================
% SOC CONSTRUCTION
% =========================================================================

function [Hunit,RhatVec,nnPairs] = build_board_nn_soc_unit(model)
% BOARD-NN SOC unit matrix; multiply by lambdaBoard = gamma_ij*t [eV].
%
% First-quantized form implemented on each physical NN bond i-j:
%
%   H_SOC,ij = i * lambdaBoard * ( S . Rhat_ij )
%   H_SOC,ji = H_SOC,ij^\dagger
%
% where S = sigma/2 and
%
%   Rhat_ij = (R_j-R_i)/|R_j-R_i| .
%
% IMPORTANT CONVENTION:
% lambdaBoard is the coefficient gamma_ij*t written on the board.
% Therefore lambdaBoard = 5 meV means that 5 meV multiplies S.Rhat.
% Since S=sigma/2, the two eigenvalues of S.Rhat are +/-1/2.
%
% The factor i gives the standard time-reversal-compatible imaginary
% spin-dependent hopping; Hermiticity is enforced by the reverse block.

L = size(model.Ri,1);
N = 2*L;
Nb = size(model.bonds,1);

sx = [0 1;1 0];
sy = [0 -1i;1i 0];
sz = [1 0;0 -1];
Sx = 0.5*sx;
Sy = 0.5*sy;
Sz = 0.5*sz;

Hunit = complex(zeros(N,N));
RhatVec = zeros(Nb,3);
nnPairs = model.bonds;

for b = 1:Nb
    i = model.bonds(b,1);
    j = model.bonds(b,2);

    dR = model.Ri(j,:) - model.Ri(i,:);
    Rhat = dR/norm(dR);
    RhatVec(b,:) = Rhat;

    SdotR = Rhat(1)*Sx + Rhat(2)*Sy + Rhat(3)*Sz;
    block = 1i*SdotR;

    idxi = (2*i-1):(2*i);
    idxj = (2*j-1):(2*j);

    Hunit(idxi,idxj) = Hunit(idxi,idxj) + block;
    Hunit(idxj,idxi) = Hunit(idxj,idxi) + block';
end

end


function caseModel = build_soc_case_model(model,soc_meV)

caseModel = model;
lambdaSOC_eV = 1e-3*soc_meV;
caseModel.lambdaSOC_meV = soc_meV;
caseModel.H_static_spin = model.H_surface_spin_noSOC + ...
    lambdaSOC_eV*model.Hsoc_unit;

% Precompute physical NN backbone hopping blocks.  BOARD-NN SOC lives on the same physical NN bonds, so the current operator
% must use the complete spin-dependent 2x2 hopping block H_ij.
Nb = size(model.bonds,1);
caseModel.bondHij = complex(zeros(2,2,Nb));

for b = 1:Nb
    i = model.bonds(b,1);
    j = model.bonds(b,2);
    idxi = (2*i-1):(2*i);
    idxj = (2*j-1):(2*j);
    caseModel.bondHij(:,:,b) = caseModel.H_static_spin(idxi,idxj);
end

end


%% =========================================================================
% PHYSICAL NN BACKBONE CURRENT (rho CAN INCLUDE SOC)
% =========================================================================

function J_eV = compute_bond_current_general(rho,model)
% Current on the PHYSICAL nearest-neighbour backbone bonds.
% The electronic state rho is evolved with BOARD-NN SOC. Because SOC lives
% on the physical NN bonds, the current operator uses the complete 2x2 bond block.  For H_ij=-t I this reduces exactly to
%   J = 2 t Im Tr(rho_ji),
% reproducing the original scalar-hopping current.

Nb = size(model.bonds,1);
J_eV = zeros(1,Nb);

for b = 1:Nb
    i = model.bonds(b,1);
    j = model.bonds(b,2);
    idxi = (2*i-1):(2*i);
    idxj = (2*j-1):(2*j);

    Hij = model.bondHij(:,:,b);
    rhoji = rho(idxj,idxi);

    J_eV(b) = -2*imag(trace(Hij*rhoji));
end

end


function J_eV = compute_bond_current_scalar_reference(rho,model)

Nb = size(model.bonds,1);
J_eV = zeros(1,Nb);

for b = 1:Nb
    i = model.bonds(b,1);
    j = model.bonds(b,2);
    tij = model.t_bonds_eV(b);

    up_i = 2*i-1;
    dn_i = 2*i;
    up_j = 2*j-1;
    dn_j = 2*j;

    J_eV(b) = 2*tij*imag(rho(up_j,up_i)+rho(dn_j,dn_i));
end

end


%% =========================================================================
% EXCHANGE
% =========================================================================

function J_site_eV = exchange_profile(J0_meV,model,p)

distance_A = vecnorm(model.Ri-p.Rstar_A,2,2);
J_site_eV = 1e-3*J0_meV*exp(-distance_A/p.r0_A);

end


function ex = precompute_exchange_matrices(J_site_eV)

L = numel(J_site_eV);
N = 2*L;

sx = [0 1;1 0];
sy = [0 -1i;1i 0];
sz = [1 0;0 -1];

Hx = complex(zeros(N,N));
Hy = complex(zeros(N,N));
Hz = complex(zeros(N,N));

for i = 1:L
    idx = (2*i-1):(2*i);
    Hx(idx,idx) = 0.5*J_site_eV(i)*sx;
    Hy(idx,idx) = 0.5*J_site_eV(i)*sy;
    Hz(idx,idx) = 0.5*J_site_eV(i)*sz;
end

ex.Hx = Hx;
ex.Hy = Hy;
ex.Hz = Hz;

end


function B_ex_T = compute_Bex_fast(rho,J_site_eV,model,p)

ud = rho(model.idxUD);
mx = 2*real(ud);
my = -2*imag(ud);
mz = real(rho(model.idxDiagUp))-real(rho(model.idxDiagDn));

J = J_site_eV(:);
weightedSpin_eV = [sum(J.*mx(:)),sum(J.*my(:)),sum(J.*mz(:))];

B_ex_T = -weightedSpin_eV/(2*p.surface_moment_eV_T);

end


%% =========================================================================
% FAST LINDBLAD BATH (SAME BATH IN ALL SOC CASES)
% =========================================================================

function bath = build_fast_lindblad_bath(H0,p)

L = size(H0,1);
N = 2*L;

[U,D] = eig(H0);
[E,idx] = sort(real(diag(D)));
U = U(:,idx);

T = kron(U,eye(2));
W = zeros(N,N);

for n = 1:L
    for m = n+1:L
        gap_eV = E(m)-E(n);
        if gap_eV <= 1e-12
            continue;
        end

        gamma_down = p.Gamma0_per_fs;
        gamma_up = p.Gamma0_per_fs*exp(-p.beta_eV_inv*gap_eV);

        for spin = 1:2
            a_n = 2*n-2+spin;
            a_m = 2*m-2+spin;
            W(a_m,a_n) = gamma_up;
            W(a_n,a_m) = gamma_down;
        end
    end
end

outRate = sum(W,1).';

bath.T = T;
bath.Tdag = T';
bath.W = W;
bath.outRate = outRate;
bath.decayMatrix = 0.5*(outRate+outRate.');

end


function drho = apply_fast_lindblad(rho,bath)

rhoE = bath.Tdag*rho*bath.T;
dE = -bath.decayMatrix.*rhoE;

pop = diag(rhoE);
gain = bath.W*pop;

N = size(rhoE,1);
diagIdx = 1:(N+1):N^2;
dE(diagIdx) = dE(diagIdx)+gain.';

drho = bath.T*dE*bath.Tdag;

end


%% =========================================================================
% OBSERVABLE INDICES / ELECTRONIC SPIN
% =========================================================================

function model = precompute_fast_indices(model)

L = size(model.Ri,1);
N = 2*L;
up = (1:2:N).';
dn = (2:2:N).';

model.up = up;
model.dn = dn;
model.idxDiagUp = sub2ind([N N],up,up);
model.idxDiagDn = sub2ind([N N],dn,dn);
model.idxUD = sub2ind([N N],up,dn);

end


function obs = density_observables_from_rho(rho,p)

N = size(rho,1);
up = 1:2:N;
down = 2:2:N;

obs.n_up = real(diag(rho(up,up))).';
obs.n_down = real(diag(rho(down,down))).';
obs.n_total = obs.n_up+obs.n_down;

rho_ud = diag(rho(up,down)).';
obs.m_x = 2*real(rho_ud);
obs.m_y = -2*imag(rho_ud);
obs.m_z = obs.n_up-obs.n_down;

obs.P_x = nan(size(obs.n_total));
obs.P_y = nan(size(obs.n_total));
obs.P_z = nan(size(obs.n_total));

valid = obs.n_total > p.density_threshold;
obs.P_x(valid) = obs.m_x(valid)./obs.n_total(valid);
obs.P_y(valid) = obs.m_y(valid)./obs.n_total(valid);
obs.P_z(valid) = obs.m_z(valid)./obs.n_total(valid);

end


%% =========================================================================
% INITIAL STATE
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
diagnostic.initialTotalMz = real(sum(diag(rho0(up,up)))-sum(diag(rho0(down,down))));

assert(diagnostic.residual < 1e-10,'Initial GS residual too large.');
assert(diagnostic.traceError < 1e-12,'Initial trace error too large.');
assert(diagnostic.hermiticityError < 1e-12,'Initial rho is not Hermitian.');
assert(abs(diagnostic.initialTotalMz) < 1e-12,'Initial molecular spin is not unpolarized.');

end


%% =========================================================================
% GEOMETRY / HOPPING
% =========================================================================

function model = build_upright_helix_model()

Nsites = 13;
sitesPerTurn = 8;
bondLength_A = 1.54;
dphi = 2*pi/sitesPerTurn;

axialRisePerSite_A = bondLength_A/sqrt(3);
transverseChord_A = bondLength_A*sqrt(2/3);
helixRadius_A = transverseChord_A/(2*sin(dphi/2));

siteIndex = (0:Nsites-1).';
phi = siteIndex*dphi;

Ri = [ ...
    helixRadius_A*(cos(phi)-1), ...
    helixRadius_A*sin(phi), ...
    1.0+axialRisePerSite_A*siteIndex];

bonds = [(1:Nsites-1).' (2:Nsites).'];
bondLengths_A = vecnorm(diff(Ri,1,1),2,2);

assert(max(abs(bondLengths_A-bondLength_A)) < 1e-12, ...
    'Helix bond lengths are not exactly 1.54 A.');

b1 = Ri(2,:)-Ri(1,:);
b2 = Ri(3,:)-Ri(2,:);
b3 = Ri(4,:)-Ri(3,:);
chiralityTripleProduct = dot(b1,cross(b2,b3));

assert(abs(chiralityTripleProduct) > 1e-8,'Helix lost 3D chirality.');

model.name = 'helical_upright_13site_8perturn_1p5turn';
model.Ri = Ri;
model.bonds = bonds;
model.t_bonds = ones(size(bonds,1),1);
model.sitesPerTurn = sitesPerTurn;
model.nTurns = (Nsites-1)/sitesPerTurn;
model.helixRadius_A = helixRadius_A;
model.pitch_A = sitesPerTurn*axialRisePerSite_A;
model.chiralityTripleProduct = chiralityTripleProduct;

end


function audit_geometry(model)

bondLengths_A = zeros(size(model.bonds,1),1);
for b = 1:size(model.bonds,1)
    bondLengths_A(b) = norm(model.Ri(model.bonds(b,2),:)-model.Ri(model.bonds(b,1),:));
end

fprintf('\nGeometry audit:\n');
fprintf('  sites              = %d\n',size(model.Ri,1));
fprintf('  bonds              = %d\n',size(model.bonds,1));
fprintf('  bond range         = %.12f ... %.12f A\n',min(bondLengths_A),max(bondLengths_A));
fprintf('  helix radius       = %.12f A\n',model.helixRadius_A);
fprintf('  turns              = %.3f\n',model.nTurns);
fprintf('  chirality triple   = %+.6e A^3\n',model.chiralityTripleProduct);

end


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


%% =========================================================================
% BIOT-SAVART
% =========================================================================

function kernel = precompute_BS_kernel_at_points(Ri,bonds,Rpoints_A,Nseg,min_distance_A)

nPoints = size(Rpoints_A,1);
Nb = size(bonds,1);
kernel = zeros(3*nPoints,Nb);

for activeBond = 1:Nb
    unit_J_eV = zeros(1,Nb);
    unit_J_eV(activeBond) = 1;
    Bunit_mT = zeros(nPoints,3);

    for ip = 1:nPoints
        Bunit_mT(ip,:) = biot_savart_at_point( ...
            Ri,bonds,unit_J_eV,Rpoints_A(ip,:),0,Nseg,min_distance_A);
    end

    kernel(:,activeBond) = reshape(Bunit_mT.',[],1);
end

end


function B_mT = apply_BS_kernel(kernel,J_eV)

Bflat_mT = kernel*J_eV(:);
nPoints = size(kernel,1)/3;
B_mT = reshape(Bflat_mT,3,nPoints).';

end


function B_mT = biot_savart_at_point(Ri,bonds,J_eV,Robs_A,excludedSite,Nseg,min_distance_A)

e_C = 1.602176634e-19;
hbar_eV_s = 6.582119569e-16;
I_A = (e_C/hbar_eV_s)*J_eV;
mu0_over_4pi = 1e-7;
B_T = [0 0 0];

for b = 1:size(bonds,1)
    if abs(I_A(b)) < 1e-30
        continue;
    end

    i = bonds(b,1);
    j = bonds(b,2);

    if excludedSite > 0 && (i==excludedSite || j==excludedSite)
        continue;
    end

    r1_A = Ri(i,:);
    r2_A = Ri(j,:);
    dl_A = (r2_A-r1_A)/Nseg;

    for s = 1:Nseg
        midpoint_A = r1_A+(s-0.5)*dl_A;
        rvec_A = Robs_A-midpoint_A;
        r_A = norm(rvec_A);

        if r_A < min_distance_A
            continue;
        end

        dl_m = 1e-10*dl_A;
        rvec_m = 1e-10*rvec_A;
        r_m = norm(rvec_m);

        dB_T = mu0_over_4pi*I_A(b)*cross(dl_m,rvec_m)/(r_m^3);
        B_T = B_T+dB_T;
    end
end

B_mT = 1e3*B_T;

end


%% =========================================================================
% OUTPUT HELPERS
% =========================================================================

function write_case_csv(E,outputFolder)

tagS = safe_number_tag(E.lambdaSOC_meV);

T = table( ...
    E.t_fs(:),E.t_fs(:)/1000, ...
    E.M(:,1),E.M(:,2),E.M(:,3), ...
    E.B_mol_mT(:,1),E.B_mol_mT(:,2),E.B_mol_mT(:,3), ...
    E.B_ex_T(:,1),E.B_ex_T(:,2),E.B_ex_T(:,3), ...
    E.B_eff_T(:,1),E.B_eff_T(:,2),E.B_eff_T(:,3), ...
    E.B_LLG_T(:,1),E.B_LLG_T(:,2),E.B_LLG_T(:,3), ...
    E.Bmol_perp_T(:),E.Bex_perp_T(:),E.Beff_perp_T(:), ...
    E.Bmol_perp_vec_T(:,1),E.Bmol_perp_vec_T(:,2),E.Bmol_perp_vec_T(:,3), ...
    E.Bex_perp_vec_T(:,1),E.Bex_perp_vec_T(:,2),E.Bex_perp_vec_T(:,3), ...
    E.Beff_perp_vec_T(:,1),E.Beff_perp_vec_T(:,2),E.Beff_perp_vec_T(:,3), ...
    E.totalElectronicSpin(:,1),E.totalElectronicSpin(:,2),E.totalElectronicSpin(:,3), ...
    E.spin_perp(:),E.spin_perp_vec(:,1),E.spin_perp_vec(:,2),E.spin_perp_vec(:,3), ...
    E.dM_mag_per_fs(:),E.dMz_per_fs(:),E.tilt_deg(:), ...
    E.traceError(:),E.hermiticityError(:),E.MnormError(:), ...
    'VariableNames',{ ...
    'time_fs','time_ps','Mx','My','Mz', ...
    'Bmol_x_mT','Bmol_y_mT','Bmol_z_mT', ...
    'Bex_x_T','Bex_y_T','Bex_z_T', ...
    'Beff_x_T','Beff_y_T','Beff_z_T', ...
    'BLLG_x_T','BLLG_y_T','BLLG_z_T', ...
    'Bmol_perp_T','Bex_perp_T','Beff_perp_T', ...
    'Bmol_perp_x_T','Bmol_perp_y_T','Bmol_perp_z_T', ...
    'Bex_perp_x_T','Bex_perp_y_T','Bex_perp_z_T', ...
    'Beff_perp_x_T','Beff_perp_y_T','Beff_perp_z_T', ...
    'electronic_mx_total','electronic_my_total','electronic_mz_total', ...
    'electronic_spin_perp','spin_perp_x','spin_perp_y','spin_perp_z', ...
    'dM_mag_per_fs','dMz_per_fs','tilt_deg', ...
    'trace_error','hermiticity_error','Mnorm_error'});

writetable(T,fullfile(outputFolder, ...
    sprintf('SOC_%s_meV_EVERY_SOLVER_STEP.csv',tagS)));

Tf = table( ...
    (1:numel(E.finalObs.m_z)).', ...
    E.finalObs.n_up(:),E.finalObs.n_down(:), ...
    E.finalObs.m_x(:),E.finalObs.m_y(:),E.finalObs.m_z(:), ...
    'VariableNames',{'site','n_up','n_down','m_x','m_y','m_z'});

writetable(Tf,fullfile(outputFolder, ...
    sprintf('SOC_%s_meV_FINAL_SITE_SPIN.csv',tagS)));

end


function save_figure(fig,outputFolder,fileName)

savefig(fig,fullfile(outputFolder,strrep(fileName,'.png','.fig')));
exportgraphics(fig,fullfile(outputFolder,fileName),'Resolution',300);

end


function tag = safe_number_tag(x)

tag = sprintf('%.9g',x);
tag = strrep(tag,'.','p');
tag = strrep(tag,'-','m');
tag = strrep(tag,'+','');
end


function downloadsDir = get_downloads_directory()

if ispc
    downloadsDir = fullfile(getenv('USERPROFILE'),'Downloads');
else
    homeDir = getenv('HOME');
    downloadsDir = fullfile(homeDir,'Downloads');
end

if ~exist(downloadsDir,'dir')
    downloadsDir = pwd;
end

end
