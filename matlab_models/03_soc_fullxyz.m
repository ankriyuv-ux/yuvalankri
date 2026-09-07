function soc_fullxyz()
% HELIX_SOC_BOARD_NN_FULLXYZ_BOTH_J0_10PS
% =========================================================================
% CLEAN 2x4 COMPARISON FOR THE NEW BOARD-NN SOC MODEL
%
% Runs BOTH missing full-vector models in one call:
%
%   A) J0 = 0 meV   : exchange OFF
%   B) J0 = 150 meV : exchange ON
%
% For EACH J0, the exact same SOC scan is used:
%   gamma*t = [0, 0.1, +5, -5] meV
%
% THE ONLY LLG CHANGE relative to the existing Z-only controls is:
%
%   OLD Z-only: B_LLG = [0,0,B_eff,z]
%   THIS RUN  : B_LLG = B_eff = [B_eff,x,B_eff,y,B_eff,z]
%
% with
%   B_eff = B_mol + B_ex.
%
% Thus the new runs are directly comparable to the two existing Z-only
% BOARD-NN data sets:
%   - J0=0   BOARD-NN Z-only
%   - J0=150 BOARD-NN Z-only
%
% BOARD-NN SOC (first quantization):
%   H_SOC,ij = i*(gamma*t)*(S . Rhat_ij)
%   S = sigma/2
%   Rhat_ij = (R_j-R_i)/|R_j-R_i|
%
% IMPORTANT:
%   * Same 13-site helix, same initial rho, same bath, same M(0)=+X.
%   * Same lambda_LLG=0.6, tau_rel=500 fs, runtime=10 ps.
%   * No direct electronic Zeeman term.
%   * ode45 Refine=1; EVERY accepted solver point is retained.
%   * No downsampling/interpolation is used in the dynamics.
%
% OUTPUT:
%   One Downloads folder containing all 8 trajectories, CSV summaries,
%   and matched comparison figures for J0=0 versus J0=150.
% =========================================================================

clc;
close all;
totalTimer = tic;

%% 1. PARAMETERS ==========================================================
p.t_hop_eV = 1.0;
p.e2_over_4pieps0_eVA = 14.4;
p.q_eff = 0.30;
p.R_contact_A = [0 0 0];
p.Rstar_A = [-1.335757 -0.292878 0.0];

p.M0 = [1 0 0];
p.M0 = p.M0/norm(p.M0);
p.gamma0_fs_T = 1.760859e-4;
p.lambda = 0.60;

% The two models to compare. This is the only model switch here.
p.J0List_meV = [0 150];
p.r0_A = 3.0;

% Same BOARD-NN SOC values used in the prior Z-only runs.
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

p.tStart_fs = 0;
p.earlyEnd_fs = 50;
p.tEnd_fs = 10000;       % 10 ps
p.chunk_fs = 5000;

p.RelTol = 1e-8;
p.AbsTol = 1e-10;
p.InitialStep_fs = 1e-4;
p.MaxStepEarly_fs = 1e-2;

p.Nseg_BS = 40;
p.min_distance_A = 0.20;
p.density_threshold = 1e-8;

% Eight independent trajectories total. Parallel affects wall time only.
p.useParallelIfAvailable = true;
p.requestedWorkers = min(4,numel(p.J0List_meV)*numel(p.lambdaSOCList_meV));

assert(isequal(p.M0,[1 0 0]),'M(0) must be +X.');
assert(abs(p.lambda-0.6)<1e-14,'lambda_LLG must remain 0.6.');
assert(abs(p.tau_rel_fs-500)<1e-12,'tau_rel must remain 500 fs.');
assert(any(p.lambdaSOCList_meV==0),'SOC=0 control is required.');

fprintf('\n============================================================\n');
fprintf('BOARD-NN SOC | FULL XYZ LLG | BOTH J0 VALUES\n');
fprintf('============================================================\n');
fprintf('J0 scan           = '); fprintf('%g ',p.J0List_meV); fprintf('meV\n');
fprintf('SOC scan          = '); fprintf('%g ',p.lambdaSOCList_meV); fprintf('meV\n');
fprintf('M(0)              = [%g %g %g]\n',p.M0);
fprintf('lambda_LLG        = %.3f\n',p.lambda);
fprintf('tau_rel           = %.1f fs\n',p.tau_rel_fs);
fprintf('physical time     = %.2f ps\n',p.tEnd_fs/1000);
fprintf('SOC model         = BOARD-NN: i*(gamma*t)*(S.Rhat_ij), S=sigma/2\n');
fprintf('LLG field         = FULL B_eff = (Bx,By,Bz)\n');
fprintf('electronic Zeeman = OFF\n');
fprintf('solver            = ode45, Refine=1, no downsampling\n');
fprintf('============================================================\n');

%% 2. GEOMETRY + STATIC HAMILTONIAN ======================================
model = build_upright_helix_model();
audit_geometry(model);

Ri = model.Ri;
bonds = model.bonds;
L = size(Ri,1);

model.t_bonds_eV = p.t_hop_eV*model.t_bonds(:);
H_hop = build_hopping_hamiltonian(L,bonds,model.t_bonds_eV);

distance_contact_A = vecnorm(Ri-p.R_contact_A,2,2);
distance_contact_A(distance_contact_A<1e-12) = 1e-12;
epsilon_site_eV = p.q_eff*p.e2_over_4pieps0_eVA./distance_contact_A;
H_surface_orbital = H_hop + diag(epsilon_site_eV);
model.H_surface_spin_noSOC = kron(H_surface_orbital,eye(2));

[model.Hsoc_unit,model.soc_Rhat,model.soc_nn_pairs] = build_board_nn_soc_unit(model);
assert(norm(model.Hsoc_unit-model.Hsoc_unit','fro')<1e-12, ...
    'H_SOC unit matrix is not Hermitian.');

rhatNorm = vecnorm(model.soc_Rhat,2,2);
fprintf('\nBOARD-NN geometry audit: links=%d | |Rhat| %.9f ... %.9f\n', ...
    size(model.soc_nn_pairs,1),min(rhatNorm),max(rhatNorm));

model = precompute_fast_indices(model);
model.BS_kernel_Rstar = precompute_BS_kernel_at_points( ...
    Ri,bonds,p.Rstar_A,p.Nseg_BS,p.min_distance_A);

%% 3. INITIAL ELECTRON STATE + COMMON BATH ===============================
[rho0,initialDiagnostic] = build_gas_ground_initial_state(H_hop);
bath = build_fast_lindblad_bath(H_surface_orbital,p);
dr0 = apply_fast_lindblad(rho0,bath);
assert(abs(trace(dr0))<1e-11,'Lindblad trace check failed.');

fprintf('Initial rho: trace err=%.3e | Herm err=%.3e | total Mz=%.3e\n', ...
    initialDiagnostic.traceError,initialDiagnostic.hermiticityError, ...
    initialDiagnostic.initialTotalMz);

%% 4. VALIDATE SOC-CAPABLE NN CURRENT AT SOC=0 ===========================
case0 = build_soc_case_model(model,0);
Jold = compute_bond_current_scalar_reference(rho0,model);
Jnew = compute_bond_current_general(rho0,case0);
assert(norm(Jold-Jnew)<1e-13, ...
    'SOC-capable current failed to reproduce SOC=0 baseline.');
fprintf('Physical NN current validation at SOC=0: PASS\n');

%% 5. OUTPUT FOLDER =======================================================
downloadsDir = get_downloads_directory();
runTag = datestr(now,'yyyymmdd_HHMMSS');
outputFolder = fullfile(downloadsDir, ...
    ['SOC_BOARD_NN_FULLXYZ_BOTH_J0_LAMBDA06_10PS_' runTag]);
if ~exist(outputFolder,'dir'), mkdir(outputFolder); end
fprintf('\nOutput folder:\n%s\n',outputFolder);

%% 6. PRECOMPUTE EXCHANGE DATA FOR EACH J0 ===============================
nJ = numel(p.J0List_meV);
nSOC = numel(p.lambdaSOCList_meV);
JsiteCell = cell(nJ,1);
exchangeCell = cell(nJ,1);
for ij = 1:nJ
    JsiteCell{ij} = exchange_profile(p.J0List_meV(ij),model,p);
    exchangeCell{ij} = precompute_exchange_matrices(JsiteCell{ij});
end

%% 7. PARALLEL AVAILABILITY ==============================================
useParallel = false;
if p.useParallelIfAvailable && license('test','Distrib_Computing_Toolbox')
    try
        pool = gcp('nocreate');
        if isempty(pool), pool = parpool('local',p.requestedWorkers); end
        useParallel = pool.NumWorkers>=2;
    catch ME
        fprintf('Parallel unavailable: %s\nRunning serially.\n',ME.message);
        useParallel = false;
    end
end

if useParallel
    fprintf('\nAll 8 full-XYZ cases will run in PARALLEL where workers allow.\n');
else
    fprintf('\nAll 8 full-XYZ cases will run SERIALLY.\n');
end

%% 8. RUN ALL J0 x SOC CASES =============================================
% runs{ij,is}: ij indexes J0, is indexes SOC.
runs = cell(nJ,nSOC);

caseJidx = zeros(nJ*nSOC,1);
caseSidx = zeros(nJ*nSOC,1);
ic = 0;
for ij=1:nJ
    for is=1:nSOC
        ic=ic+1;
        caseJidx(ic)=ij;
        caseSidx(ic)=is;
    end
end

flatRuns = cell(numel(caseJidx),1);
if useParallel
    parfor icase = 1:numel(caseJidx)
        ij = caseJidx(icase);
        is = caseSidx(icase);
        flatRuns{icase} = run_one_full_case( ...
            rho0,p.J0List_meV(ij),p.lambdaSOCList_meV(is), ...
            JsiteCell{ij},exchangeCell{ij},model,bath,p);
    end
else
    for icase = 1:numel(caseJidx)
        ij = caseJidx(icase);
        is = caseSidx(icase);
        flatRuns{icase} = run_one_full_case( ...
            rho0,p.J0List_meV(ij),p.lambdaSOCList_meV(is), ...
            JsiteCell{ij},exchangeCell{ij},model,bath,p);
    end
end

for icase=1:numel(caseJidx)
    runs{caseJidx(icase),caseSidx(icase)} = flatRuns{icase};
end

%% 9. SAVE CASES + BUILD ONE CLEAN SUMMARY TABLE =========================
rows = nJ*nSOC;
J0_col = zeros(rows,1);
soc_col = zeros(rows,1);
Mz_end = zeros(rows,1);
My_end = zeros(rows,1);
Mx_end = zeros(rows,1);
maxAbsMz = zeros(rows,1);
maxTilt_deg = zeros(rows,1);
maxBmolPerp_mT = zeros(rows,1);
maxBexPerp_mT = zeros(rows,1);
maxBeffPerp_mT = zeros(rows,1);
maxSpinPerp = zeros(rows,1);
maxTorque_per_fs = zeros(rows,1);
maxAbs_dMz_per_fs = zeros(rows,1);
maxAbsBmolZ_mT = zeros(rows,1);
maxAbsBexZ_mT = zeros(rows,1);
maxAbsBeffZ_mT = zeros(rows,1);
signedBeffZArea_Tfs = zeros(rows,1);
absBeffZArea_Tfs = zeros(rows,1);
cancelBeffZ_end = zeros(rows,1);
nSolverPoints = zeros(rows,1);
runtime_min = zeros(rows,1);

ir=0;
for ij=1:nJ
    for is=1:nSOC
        ir=ir+1;
        E=runs{ij,is};

        % Z-area remains a useful cancellation diagnostic, but unlike the
        % Z-only runs it is NOT an exact predictor of Mz under full XYZ.
        BeffZ_T = E.B_eff_T(:,3);
        Az_signed = cumtrapz(E.t_fs(:),BeffZ_T);
        Az_abs = cumtrapz(E.t_fs(:),abs(BeffZ_T));
        E.cumBeffZ_signed_Tfs = Az_signed;
        E.cumBeffZ_abs_Tfs = Az_abs;
        E.BeffZ_cancellationRatio = abs(Az_signed)./max(Az_abs,eps);
        E.BeffZ_cancellationRatio(1)=NaN;
        runs{ij,is}=E;

        J0_col(ir)=E.J0_meV;
        soc_col(ir)=E.lambdaSOC_meV;
        Mx_end(ir)=E.M(end,1);
        My_end(ir)=E.M(end,2);
        Mz_end(ir)=E.M(end,3);
        maxAbsMz(ir)=max(abs(E.M(:,3)));
        maxTilt_deg(ir)=max(E.tilt_deg);
        maxBmolPerp_mT(ir)=1e3*max(E.Bmol_perp_T);
        maxBexPerp_mT(ir)=1e3*max(E.Bex_perp_T);
        maxBeffPerp_mT(ir)=1e3*max(E.Beff_perp_T);
        maxSpinPerp(ir)=max(E.spin_perp);
        maxTorque_per_fs(ir)=max(E.dM_mag_per_fs);
        maxAbs_dMz_per_fs(ir)=max(abs(E.dMz_per_fs));
        maxAbsBmolZ_mT(ir)=max(abs(E.B_mol_mT(:,3)));
        maxAbsBexZ_mT(ir)=1e3*max(abs(E.B_ex_T(:,3)));
        maxAbsBeffZ_mT(ir)=1e3*max(abs(E.B_eff_T(:,3)));
        signedBeffZArea_Tfs(ir)=Az_signed(end);
        absBeffZArea_Tfs(ir)=Az_abs(end);
        cancelBeffZ_end(ir)=abs(Az_signed(end))/max(Az_abs(end),eps);
        nSolverPoints(ir)=numel(E.t_fs);
        runtime_min(ir)=E.runtime_min;

        tagJ=safe_number_tag(E.J0_meV);
        tagS=safe_number_tag(E.lambdaSOC_meV);
        save(fullfile(outputFolder,sprintf('J0_%s_SOC_%s_FULLXYZ.mat',tagJ,tagS)), ...
            'E','-v7.3');
        write_case_csv(E,outputFolder);
    end
end

summaryTable = table(J0_col,soc_col,Mx_end,My_end,Mz_end,maxAbsMz,maxTilt_deg, ...
    maxBmolPerp_mT,maxBexPerp_mT,maxBeffPerp_mT,maxSpinPerp, ...
    maxTorque_per_fs,maxAbs_dMz_per_fs,maxAbsBmolZ_mT,maxAbsBexZ_mT, ...
    maxAbsBeffZ_mT,signedBeffZArea_Tfs,absBeffZArea_Tfs,cancelBeffZ_end, ...
    nSolverPoints,runtime_min, ...
    'VariableNames',{'J0_meV','lambda_SOC_meV','Mx_end','My_end','Mz_end', ...
    'max_abs_Mz','max_tilt_deg','max_Bmol_perp_mT','max_Bex_perp_mT', ...
    'max_Beff_perp_mT','max_electronic_spin_perp','max_dM_dt_per_fs', ...
    'max_abs_dMz_dt_per_fs','max_abs_Bmol_z_mT','max_abs_Bex_z_mT', ...
    'max_abs_Beff_z_mT','signed_Beff_z_area_Tfs','abs_Beff_z_area_Tfs', ...
    'Beff_z_cancellation_ratio_end','n_solver_points','runtime_min'});

writetable(summaryTable,fullfile(outputFolder,'BOARD_NN_FULLXYZ_BOTH_J0_SUMMARY.csv'));
save(fullfile(outputFolder,'COMPLETE_BOARD_NN_FULLXYZ_BOTH_J0.mat'), ...
    'runs','summaryTable','p','model','initialDiagnostic','-v7.3');

fprintf('\n============================================================\n');
fprintf('FULL XYZ SUMMARY: BOARD-NN SOC, BOTH J0 VALUES\n');
fprintf('============================================================\n');
disp(summaryTable);

%% 10. FIGURES ============================================================
labels = arrayfun(@(x)sprintf('SOC = %+g meV',x),p.lambdaSOCList_meV, ...
    'UniformOutput',false);

% FIGURE 0: geometry ------------------------------------------------------
fig0=figure('Color','w','Position',[60 60 1200 850]);
plot3(model.Ri(:,1),model.Ri(:,2),model.Ri(:,3),'-o','LineWidth',1.2); hold on;
for ii=1:size(model.soc_nn_pairs,1)
    i=model.soc_nn_pairs(ii,1); j=model.soc_nn_pairs(ii,2);
    rmid=0.5*(model.Ri(i,:)+model.Ri(j,:)); rh=model.soc_Rhat(ii,:);
    quiver3(rmid(1),rmid(2),rmid(3),rh(1),rh(2),rh(3),0.55,'LineWidth',1.0);
end
grid on; axis equal; xlabel('x [A]'); ylabel('y [A]'); zlabel('z [A]');
title('BOARD-NN SOC geometry: Rhat_{ij} on physical NN bonds'); view(35,25);
save_figure(fig0,outputFolder,'FIG00_BOARD_NN_SOC_GEOMETRY.png');

% FIGURE 1: Mz(t), direct J0 comparison ----------------------------------
fig1=figure('Color','w','Position',[80 70 1500 900]);
tiledlayout(nJ,1,'TileSpacing','compact','Padding','compact');
for ij=1:nJ
    nexttile; hold on;
    for is=1:nSOC
        E=runs{ij,is};
        plot(E.t_fs/1000,E.M(:,3),'LineWidth',1.3,'DisplayName',labels{is});
    end
    grid on; ylabel('M_z'); title(sprintf('J_0 = %g meV | FULL B_{eff}',p.J0List_meV(ij)));
    legend('Location','best'); xlim([0 p.tEnd_fs/1000]);
    if ij==nJ, xlabel('time [ps]'); end
end
sgtitle('BOARD-NN SOC: out-of-plane moment response under full XYZ LLG');
save_figure(fig1,outputFolder,'FIG01_MZ_FULLXYZ_BOTH_J0.png');

% FIGURE 2: all M components for +/-5, separated by J0 -------------------
fig2=figure('Color','w','Position',[100 50 1600 1050]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
for ij=1:nJ
    ip=find(abs(p.lambdaSOCList_meV-5)<1e-12,1);
    im=find(abs(p.lambdaSOCList_meV+5)<1e-12,1);
    Ep=runs{ij,ip}; Em=runs{ij,im};
    nexttile((ij-1)*3+1); hold on;
    plot(Ep.t_fs/1000,Ep.M(:,1),'LineWidth',1.2,'DisplayName','+5');
    plot(Em.t_fs/1000,Em.M(:,1),'LineWidth',1.2,'DisplayName','-5');
    grid on; ylabel('M_x'); title(sprintf('J_0=%g meV',p.J0List_meV(ij))); legend('Location','best');
    nexttile((ij-1)*3+2); hold on;
    plot(Ep.t_fs/1000,Ep.M(:,2),'LineWidth',1.2,'DisplayName','+5');
    plot(Em.t_fs/1000,Em.M(:,2),'LineWidth',1.2,'DisplayName','-5');
    grid on; ylabel('M_y'); title('M_y');
    nexttile((ij-1)*3+3); hold on;
    plot(Ep.t_fs/1000,Ep.M(:,3),'LineWidth',1.2,'DisplayName','+5');
    plot(Em.t_fs/1000,Em.M(:,3),'LineWidth',1.2,'DisplayName','-5');
    grid on; ylabel('M_z'); title('M_z');
    if ij==nJ
        for kk=1:3, nexttile((ij-1)*3+kk); xlabel('time [ps]'); end
    end
end
sgtitle('Moment components under full XYZ field | BOARD-NN SOC = +/-5 meV');
save_figure(fig2,outputFolder,'FIG02_M_COMPONENTS_PLUS_MINUS5_BOTH_J0.png');

% FIGURE 3: full Beff components for +5 meV, both J0 ---------------------
fig3=figure('Color','w','Position',[120 60 1550 950]);
tiledlayout(nJ,1,'TileSpacing','compact','Padding','compact');
ip=find(abs(p.lambdaSOCList_meV-5)<1e-12,1);
for ij=1:nJ
    E=runs{ij,ip}; nexttile; hold on;
    plot(E.t_fs/1000,E.B_eff_T(:,1),'LineWidth',1.1,'DisplayName','B_{eff,x}');
    plot(E.t_fs/1000,E.B_eff_T(:,2),'LineWidth',1.1,'DisplayName','B_{eff,y}');
    plot(E.t_fs/1000,E.B_eff_T(:,3),'LineWidth',1.1,'DisplayName','B_{eff,z}');
    grid on; ylabel('B_{eff} [T]'); title(sprintf('J_0=%g meV | SOC=+5 meV',p.J0List_meV(ij)));
    legend('Location','best'); if ij==nJ, xlabel('time [ps]'); end
end
sgtitle('Full effective-field components actually entering LLG');
save_figure(fig3,outputFolder,'FIG03_BEFF_COMPONENTS_PLUS5_BOTH_J0.png');

% FIGURE 4: exchange field components at J0=150, +5 meV ------------------
fig4=figure('Color','w','Position',[140 80 1500 850]);
E=runs{find(abs(p.J0List_meV-150)<1e-12,1),ip}; hold on;
plot(E.t_fs/1000,E.B_ex_T(:,1),'LineWidth',1.15,'DisplayName','B_{ex,x}');
plot(E.t_fs/1000,E.B_ex_T(:,2),'LineWidth',1.15,'DisplayName','B_{ex,y}');
plot(E.t_fs/1000,E.B_ex_T(:,3),'LineWidth',1.15,'DisplayName','B_{ex,z}');
grid on; xlabel('time [ps]'); ylabel('B_{ex} [T]'); legend('Location','best');
title('Exchange-field components | J_0=150 meV | BOARD-NN SOC=+5 meV | full XYZ');
save_figure(fig4,outputFolder,'FIG04_BEX_COMPONENTS_J0150_PLUS5.png');

% FIGURE 5: transverse effective field -----------------------------------
fig5=figure('Color','w','Position',[160 70 1500 900]);
tiledlayout(nJ,1,'TileSpacing','compact','Padding','compact');
for ij=1:nJ
    nexttile; hold on;
    for is=1:nSOC
        E=runs{ij,is};
        plot(E.t_fs/1000,1e3*E.Beff_perp_T,'LineWidth',1.25,'DisplayName',labels{is});
    end
    grid on; ylabel('|B_{eff,\perp}| [mT]'); title(sprintf('J_0=%g meV',p.J0List_meV(ij)));
    legend('Location','best'); if ij==nJ, xlabel('time [ps]'); end
end
sgtitle('Transverse effective field relative to instantaneous M');
save_figure(fig5,outputFolder,'FIG05_BEFF_PERP_BOTH_J0.png');

% FIGURE 6: transverse electronic spin -----------------------------------
fig6=figure('Color','w','Position',[180 70 1500 900]);
tiledlayout(nJ,1,'TileSpacing','compact','Padding','compact');
for ij=1:nJ
    nexttile; hold on;
    for is=1:nSOC
        E=runs{ij,is};
        plot(E.t_fs/1000,E.spin_perp,'LineWidth',1.25,'DisplayName',labels{is});
    end
    grid on; ylabel('|s_{e,\perp}|'); title(sprintf('J_0=%g meV',p.J0List_meV(ij)));
    legend('Location','best'); if ij==nJ, xlabel('time [ps]'); end
end
sgtitle('Electronic transverse-spin response');
save_figure(fig6,outputFolder,'FIG06_ELECTRON_SPIN_PERP_BOTH_J0.png');

% FIGURE 7: actual full-vector LLG torque --------------------------------
fig7=figure('Color','w','Position',[200 70 1500 900]);
tiledlayout(nJ,1,'TileSpacing','compact','Padding','compact');
for ij=1:nJ
    nexttile; hold on;
    for is=1:nSOC
        E=runs{ij,is};
        plot(E.t_fs/1000,E.dM_mag_per_fs,'LineWidth',1.25,'DisplayName',labels{is});
    end
    grid on; ylabel('|dM/dt| [fs^{-1}]'); title(sprintf('J_0=%g meV',p.J0List_meV(ij)));
    legend('Location','best'); if ij==nJ, xlabel('time [ps]'); end
end
sgtitle('Actual LLG rotation rate under full B_{eff}');
save_figure(fig7,outputFolder,'FIG07_LLG_TORQUE_BOTH_J0.png');

% FIGURE 8: compact numerical summary vs SOC -----------------------------
fig8=figure('Color','w','Position',[220 70 1500 950]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
styles={'-o','-s'};
for metric=1:4
    nexttile; hold on;
    for ij=1:nJ
        vals=zeros(1,nSOC);
        for is=1:nSOC
            E=runs{ij,is};
            switch metric
                case 1, vals(is)=max(abs(E.M(:,3)));
                case 2, vals(is)=max(E.tilt_deg);
                case 3, vals(is)=1e3*max(E.Beff_perp_T);
                case 4, vals(is)=1e3*max(E.Bex_perp_T);
            end
        end
        [xs,ord]=sort(p.lambdaSOCList_meV);
        plot(xs,vals(ord),styles{ij},'LineWidth',1.35, ...
            'DisplayName',sprintf('J_0=%g meV',p.J0List_meV(ij)));
    end
    grid on; xlabel('\gamma t / SOC [meV]'); legend('Location','best');
    switch metric
        case 1, ylabel('max |M_z|'); title('Out-of-plane response');
        case 2, ylabel('max tilt from +X [deg]'); title('Moment reorientation');
        case 3, ylabel('max |B_{eff,\perp}| [mT]'); title('Total transverse field');
        case 4, ylabel('max |B_{ex,\perp}| [mT]'); title('Exchange transverse field');
    end
end
sgtitle('BOARD-NN full-XYZ summary: isolating the effect of exchange');
save_figure(fig8,outputFolder,'FIG08_FULLXYZ_SUMMARY_BOTH_J0.png');

% FIGURE 9: sign control in final Mz -------------------------------------
fig9=figure('Color','w','Position',[240 100 1300 760]); hold on;
for ij=1:nJ
    vals=zeros(1,nSOC);
    for is=1:nSOC, vals(is)=runs{ij,is}.M(end,3); end
    i0=find(abs(p.lambdaSOCList_meV)<1e-15,1);
    vals=vals-vals(i0);
    [xs,ord]=sort(p.lambdaSOCList_meV);
    plot(xs,vals(ord),styles{ij},'LineWidth',1.4, ...
        'DisplayName',sprintf('J_0=%g meV',p.J0List_meV(ij)));
end
yline(0,'--'); grid on; xlabel('\gamma t / SOC [meV]'); ylabel('\Delta M_z(end) vs SOC=0');
title('SOC sign control under full XYZ dynamics'); legend('Location','best');
save_figure(fig9,outputFolder,'FIG09_SOC_SIGN_CONTROL_BOTH_J0.png');

% FIGURE 10: z-channel cancellation remains a diagnostic -----------------
fig10=figure('Color','w','Position',[260 80 1500 900]);
tiledlayout(nJ,1,'TileSpacing','compact','Padding','compact');
for ij=1:nJ
    nexttile; hold on;
    for is=1:nSOC
        E=runs{ij,is};
        plot(E.t_fs/1000,E.BeffZ_cancellationRatio,'LineWidth',1.2,'DisplayName',labels{is});
    end
    grid on; ylabel('|\int B_zdt| / \int|B_z|dt');
    title(sprintf('J_0=%g meV | z-channel diagnostic only',p.J0List_meV(ij)));
    legend('Location','best'); if ij==nJ, xlabel('time [ps]'); end
end
sgtitle('Temporal cancellation of B_{eff,z} (diagnostic; full XYZ drives LLG)');
save_figure(fig10,outputFolder,'FIG10_BEFF_Z_CANCELLATION_DIAGNOSTIC_BOTH_J0.png');


%% 10B. FOCUSED COMPARISON: EXCHANGE + FULL XYZ, NO SOC VS SOC ==========
% Clean comparison requested for the presentation:
%   SAME J0 = 150 meV
%   SAME M(0) = +X
%   SAME Full-XYZ LLG
%   SAME probe R*
%   SAME bath / damping / runtime
%   ONLY SOC changes: 0, +5, -5 meV
%
% Therefore:
%   SOC = 0    -> exchange + Full XYZ WITHOUT SOC
%   SOC = +5   -> exchange + Full XYZ WITH positive SOC
%   SOC = -5   -> exchange + Full XYZ WITH negative SOC

ij150 = find(abs(p.J0List_meV-150)<1e-12,1);
iSOC0 = find(abs(p.lambdaSOCList_meV-0)<1e-12,1);
iSOCp = find(abs(p.lambdaSOCList_meV-5)<1e-12,1);
iSOCm = find(abs(p.lambdaSOCList_meV+5)<1e-12,1);

assert(~isempty(ij150) && ~isempty(iSOC0) && ~isempty(iSOCp) && ~isempty(iSOCm), ...
    'Focused comparison requires J0=150 and SOC = 0,+5,-5 meV.');

E0 = runs{ij150,iSOC0};
Ep = runs{ij150,iSOCp};
Em = runs{ij150,iSOCm};

% FIGURE 11A: field comparison -------------------------------------------
fig11a = figure('Color','w','Position',[260 40 1650 1050]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile; hold on;
plot(E0.t_fs/1000,E0.B_eff_T(:,1),'LineWidth',1.15,'DisplayName','B_{eff,x}');
plot(E0.t_fs/1000,E0.B_eff_T(:,2),'LineWidth',1.15,'DisplayName','B_{eff,y}');
plot(E0.t_fs/1000,E0.B_eff_T(:,3),'LineWidth',1.15,'DisplayName','B_{eff,z}');
grid on; xlabel('time [ps]'); ylabel('B_{eff} [T]');
title('NO SOC: J_0=150 meV, Full XYZ');
legend('Location','best');

nexttile; hold on;
plot(Ep.t_fs/1000,Ep.B_eff_T(:,1),'LineWidth',1.15,'DisplayName','B_{eff,x}');
plot(Ep.t_fs/1000,Ep.B_eff_T(:,2),'LineWidth',1.15,'DisplayName','B_{eff,y}');
plot(Ep.t_fs/1000,Ep.B_eff_T(:,3),'LineWidth',1.15,'DisplayName','B_{eff,z}');
grid on; xlabel('time [ps]'); ylabel('B_{eff} [T]');
title('SOC = +5 meV: J_0=150 meV, Full XYZ');
legend('Location','best');

nexttile; hold on;
plot(Em.t_fs/1000,Em.B_eff_T(:,1),'LineWidth',1.15,'DisplayName','B_{eff,x}');
plot(Em.t_fs/1000,Em.B_eff_T(:,2),'LineWidth',1.15,'DisplayName','B_{eff,y}');
plot(Em.t_fs/1000,Em.B_eff_T(:,3),'LineWidth',1.15,'DisplayName','B_{eff,z}');
grid on; xlabel('time [ps]'); ylabel('B_{eff} [T]');
title('SOC = -5 meV: J_0=150 meV, Full XYZ');
legend('Location','best');

nexttile; hold on;
plot(E0.t_fs/1000,1e3*E0.Beff_perp_T,'LineWidth',1.35,'DisplayName','SOC = 0');
plot(Ep.t_fs/1000,1e3*Ep.Beff_perp_T,'LineWidth',1.35,'DisplayName','SOC = +5 meV');
plot(Em.t_fs/1000,1e3*Em.Beff_perp_T,'LineWidth',1.35,'DisplayName','SOC = -5 meV');
grid on; xlabel('time [ps]'); ylabel('|B_{eff,\perp}| [mT]');
title('Direct transverse-field comparison');
legend('Location','best');

sgtitle('Exchange + Full XYZ: field comparison without SOC vs with SOC');
save_figure(fig11a,outputFolder,'FIG11A_J0150_FIELD_NO_SOC_VS_SOC.png');

% FIGURE 11B: Mz comparison ----------------------------------------------
fig11b = figure('Color','w','Position',[300 100 1350 780]); hold on;
plot(E0.t_fs/1000,E0.M(:,3),'LineWidth',1.5,'DisplayName','SOC = 0');
plot(Ep.t_fs/1000,Ep.M(:,3),'LineWidth',1.5,'DisplayName','SOC = +5 meV');
plot(Em.t_fs/1000,Em.M(:,3),'LineWidth',1.5,'DisplayName','SOC = -5 meV');
yline(0,'--');
grid on;
xlabel('time [ps]');
ylabel('M_z');
title('J_0=150 meV | Full XYZ | M_z: no SOC vs SOC');
legend('Location','best');
xlim([0 p.tEnd_fs/1000]);
save_figure(fig11b,outputFolder,'FIG11B_J0150_MZ_NO_SOC_VS_SOC.png');

% Small numeric printout for presentation --------------------------------
fprintf('\n============================================================\n');
fprintf('FOCUSED J0=150 FULL-XYZ COMPARISON: NO SOC VS SOC\n');
fprintf('============================================================\n');
fprintf('SOC=0   : Mz(end)=%+.6e | max|Mz|=%.6e | max|Beff_perp|=%.6f mT\n', ...
    E0.M(end,3),max(abs(E0.M(:,3))),1e3*max(E0.Beff_perp_T));
fprintf('SOC=+5  : Mz(end)=%+.6e | max|Mz|=%.6e | max|Beff_perp|=%.6f mT\n', ...
    Ep.M(end,3),max(abs(Ep.M(:,3))),1e3*max(Ep.Beff_perp_T));
fprintf('SOC=-5  : Mz(end)=%+.6e | max|Mz|=%.6e | max|Beff_perp|=%.6f mT\n', ...
    Em.M(end,3),max(abs(Em.M(:,3))),1e3*max(Em.Beff_perp_T));
fprintf('============================================================\n');

Tfocus = table( ...
    [0;5;-5], ...
    [E0.M(end,3);Ep.M(end,3);Em.M(end,3)], ...
    [max(abs(E0.M(:,3)));max(abs(Ep.M(:,3)));max(abs(Em.M(:,3)))], ...
    1e3*[max(E0.Beff_perp_T);max(Ep.Beff_perp_T);max(Em.Beff_perp_T)], ...
    'VariableNames',{'SOC_meV','Mz_end','max_abs_Mz','max_Beff_perp_mT'});

writetable(Tfocus,fullfile(outputFolder, ...
    'FOCUSED_J0150_NO_SOC_VS_SOC_SUMMARY.csv'));


%% 11. OPTIONAL: compare against existing Z-only MAT files ===============
% This section does NOT alter the new runs. If the previous Z-only files
% exist anywhere under Downloads, one additional figure is created.
try
    f0=dir(fullfile(downloadsDir,'**','COMPLETE_J0ZERO_BOARD_NN_SOC_BZONLY.mat'));
    f1=dir(fullfile(downloadsDir,'**','COMPLETE_BOARD_NN_EXCHANGE_BZONLY.mat'));
    if ~isempty(f0) && ~isempty(f1)
        [~,k0]=max([f0.datenum]); [~,k1]=max([f1.datenum]);
        Z0=load(fullfile(f0(k0).folder,f0(k0).name),'runs');
        Z1=load(fullfile(f1(k1).folder,f1(k1).name),'runs');
        zRuns={Z0.runs,Z1.runs};

        fig11=figure('Color','w','Position',[280 60 1550 950]);
        tiledlayout(nJ,1,'TileSpacing','compact','Padding','compact');
        for ij=1:nJ
            nexttile; hold on;
            for is=1:nSOC
                Ef=runs{ij,is}; Ez=zRuns{ij}{is};
                plot(Ef.t_fs/1000,Ef.M(:,3),'LineWidth',1.25, ...
                    'DisplayName',[labels{is} ' full XYZ']);
                plot(Ez.t_fs/1000,Ez.M(:,3),'--','LineWidth',1.0, ...
                    'DisplayName',[labels{is} ' Z-only']);
            end
            grid on; ylabel('M_z'); title(sprintf('J_0=%g meV',p.J0List_meV(ij)));
            legend('Location','best'); if ij==nJ, xlabel('time [ps]'); end
        end
        sgtitle('Direct comparison: existing Z-only versus new full-XYZ BOARD-NN runs');
        save_figure(fig11,outputFolder,'FIG11_DIRECT_ZONLY_VS_FULLXYZ.png');
        fprintf('\nPrevious Z-only MAT files found: direct comparison figure created.\n');
    else
        fprintf('\nPrevious Z-only MAT files not found automatically; skipping direct comparison figure.\n');
    end
catch ME
    fprintf('\nOptional Z-only comparison skipped: %s\n',ME.message);
end

%% 12. MODEL NOTE =========================================================
fid=fopen(fullfile(outputFolder,'MODEL_NOTE_FULLXYZ_BOTH_J0.txt'),'w');
if fid>=0
    fprintf(fid,'BOARD-NN SOC FULL-XYZ comparison, J0=0 and J0=150 meV.\n');
    fprintf(fid,'SOC values: 0, 0.1, +5, -5 meV.\n');
    fprintf(fid,'H_SOC,ij = i*(gamma*t)*(S.Rhat_ij), S=sigma/2.\n');
    fprintf(fid,'LLG receives full B_eff=(Bx,By,Bz), where B_eff=B_mol+B_ex.\n');
    fprintf(fid,'For J0=0, B_ex is identically zero.\n');
    fprintf(fid,'For J0=150, full self-consistent H_MS[M] and B_ex are active.\n');
    fprintf(fid,'No direct electronic Zeeman. ode45 Refine=1; no downsampling.\n');
    fprintf(fid,'B_eff,z cancellation metrics are diagnostic only under full XYZ.\n');
    fclose(fid);
end

fprintf('\n============================================================\n');
fprintf('BOARD-NN FULL-XYZ BOTH-J0 RUN COMPLETE\n');
fprintf('Output folder:\n%s\n',outputFolder);
fprintf('Total wall time = %.2f min\n',toc(totalTimer)/60);
fprintf('============================================================\n');

end


%% =========================================================================
% ONE FULL-XYZ CASE
% =========================================================================
function E = run_one_full_case(rho0,J0_meV,soc_meV,J_site_eV,exchange,model,bath,p)
caseTimer=tic;
fprintf('\n------------------------------------------------------------\n');
fprintf('START J0=%g meV | SOC=%g meV | FULL XYZ\n',J0_meV,soc_meV);
fprintf('------------------------------------------------------------\n');

caseModel=build_soc_case_model(model,soc_meV);
Nspin=size(caseModel.H_static_spin,1);
Nrho=Nspin^2;
rhs=@(~,y) coupled_soc_exchange_rhs_full(y,J_site_eV,exchange,caseModel,bath,p);
y0=[rho0(:);p.M0(:)];

% All accepted ode45 points are accumulated.
tObs=[]; MObs=[]; BmolObs=[]; BexObs=[]; BeffObs=[]; spinObs=[];
BmolPerpObs=[]; BexPerpObs=[]; BeffPerpObs=[]; spinPerpObs=[];
dMmagObs=[]; dMzObs=[]; tiltObs=[]; traceObs=[]; hermObs=[]; MnormObs=[];

optsEarly=odeset('InitialStep',p.InitialStep_fs,'MaxStep',p.MaxStepEarly_fs, ...
    'RelTol',p.RelTol,'AbsTol',p.AbsTol,'Refine',1,'Stats','off');
[tSeg,YSeg]=ode45(rhs,[p.tStart_fs p.earlyEnd_fs],y0,optsEarly);
[tObs,MObs,BmolObs,BexObs,BeffObs,spinObs, ...
 BmolPerpObs,BexPerpObs,BeffPerpObs,spinPerpObs, ...
 dMmagObs,dMzObs,tiltObs,traceObs,hermObs,MnormObs,lastState] = ...
    append_full_segment(tSeg,YSeg,false,tObs,MObs,BmolObs,BexObs,BeffObs,spinObs, ...
    BmolPerpObs,BexPerpObs,BeffPerpObs,spinPerpObs,dMmagObs,dMzObs,tiltObs, ...
    traceObs,hermObs,MnormObs,J_site_eV,caseModel,p);

optsLate=odeset('RelTol',p.RelTol,'AbsTol',p.AbsTol,'Refine',1,'Stats','off');
currentTime=p.earlyEnd_fs;
while currentTime<p.tEnd_fs-1e-12
    nextEnd=min(currentTime+p.chunk_fs,p.tEnd_fs);
    segTimer=tic;
    [tSeg,YSeg]=ode45(rhs,[currentTime nextEnd],lastState,optsLate);
    [tObs,MObs,BmolObs,BexObs,BeffObs,spinObs, ...
     BmolPerpObs,BexPerpObs,BeffPerpObs,spinPerpObs, ...
     dMmagObs,dMzObs,tiltObs,traceObs,hermObs,MnormObs,lastState] = ...
        append_full_segment(tSeg,YSeg,true,tObs,MObs,BmolObs,BexObs,BeffObs,spinObs, ...
        BmolPerpObs,BexPerpObs,BeffPerpObs,spinPerpObs,dMmagObs,dMzObs,tiltObs, ...
        traceObs,hermObs,MnormObs,J_site_eV,caseModel,p);
    currentTime=nextEnd;
    dt=diff(tSeg); if isempty(dt), dtMed=NaN; else, dtMed=median(dt); end
    fprintf('  J0=%g SOC=%g: reached %5.2f ps | points=%7d | median dt=%g fs | %.2f min\n', ...
        J0_meV,soc_meV,currentTime/1000,numel(tSeg),dtMed,toc(segTimer)/60);
end

rhoFinal=reshape(lastState(1:Nrho),Nspin,Nspin);
rhoFinal=0.5*(rhoFinal+rhoFinal');

E.J0_meV=J0_meV;
E.lambdaSOC_meV=soc_meV;
E.t_fs=tObs;
E.M=MObs;
E.B_mol_mT=BmolObs;
E.B_ex_T=BexObs;
E.B_eff_T=BeffObs;
E.B_LLG_T=BeffObs;  % FULL XYZ is the actual LLG field in this run.
E.totalElectronicSpin=spinObs;
E.Bmol_perp_T=BmolPerpObs;
E.Bex_perp_T=BexPerpObs;
E.Beff_perp_T=BeffPerpObs;
E.spin_perp=spinPerpObs;
E.Bmol_perp_vec_T=(1e-3*BmolObs)-sum((1e-3*BmolObs).*MObs,2).*MObs;
E.Bex_perp_vec_T=BexObs-sum(BexObs.*MObs,2).*MObs;
E.Beff_perp_vec_T=BeffObs-sum(BeffObs.*MObs,2).*MObs;
E.spin_perp_vec=spinObs-sum(spinObs.*MObs,2).*MObs;
E.dM_mag_per_fs=dMmagObs;
E.dMz_per_fs=dMzObs;
E.tilt_deg=tiltObs;
E.traceError=traceObs;
E.hermiticityError=hermObs;
E.MnormError=MnormObs;
E.finalRho=rhoFinal;
E.finalObs=density_observables_from_rho(rhoFinal,p);
E.runtime_min=toc(caseTimer)/60;

fprintf('COMPLETE J0=%g SOC=%g | %.2f min | final M=[%.6g %.6g %.6g] | max Beff_perp=%.6e T\n', ...
    J0_meV,soc_meV,E.runtime_min,E.M(end,1),E.M(end,2),E.M(end,3),max(E.Beff_perp_T));
end


%% =========================================================================
% COUPLED RHS: FULL B_eff DRIVES LLG
% =========================================================================
function dy = coupled_soc_exchange_rhs_full(y,J_site_eV,exchange,model,bath,p)
Nspin=size(model.H_static_spin,1);
Nrho=Nspin^2;
rho=reshape(y(1:Nrho),Nspin,Nspin);
rho=0.5*(rho+rho');
M=real(y(Nrho+1:Nrho+3)); M=M(:).'; M=M/max(norm(M),1e-15);

J_eV=compute_bond_current_general(rho,model);
B_mol_mT=apply_BS_kernel(model.BS_kernel_Rstar,J_eV);
B_mol_T=1e-3*B_mol_mT;

H_exchange=M(1)*exchange.Hx+M(2)*exchange.Hy+M(3)*exchange.Hz;
H_total=model.H_static_spin+H_exchange;  % no direct electronic Zeeman

drho_H=-1i/p.hbar_eV_fs*(H_total*rho-rho*H_total);
drho_D=apply_fast_lindblad(rho,bath);
drho=drho_H+drho_D;

B_ex_T=compute_Bex_fast(rho,J_site_eV,model,p);
B_eff_T=B_mol_T+B_ex_T;

% THIS is the key comparison change: ALL THREE components enter LLG.
B_LLG_T=B_eff_T;
precession=cross(M,B_LLG_T);
damping=cross(M,cross(M,B_LLG_T));
dM=p.gamma0_fs_T*(precession-p.lambda*damping);

dy=[drho(:);dM(:)];
end


%% =========================================================================
% APPEND EVERY ACCEPTED SOLVER POINT: FULL XYZ
% =========================================================================
function [tObs,MObs,BmolObs,BexObs,BeffObs,spinObs, ...
          BmolPerpObs,BexPerpObs,BeffPerpObs,spinPerpObs, ...
          dMmagObs,dMzObs,tiltObs,traceObs,hermObs,MnormObs,lastState] = ...
    append_full_segment(tSeg,YSeg,dropFirst,tObs,MObs,BmolObs,BexObs,BeffObs,spinObs, ...
    BmolPerpObs,BexPerpObs,BeffPerpObs,spinPerpObs,dMmagObs,dMzObs,tiltObs, ...
    traceObs,hermObs,MnormObs,J_site_eV,model,p)

if dropFirst
    tSeg=tSeg(2:end); YSeg=YSeg(2:end,:);
end
n=numel(tSeg);
Nspin=size(model.H_static_spin,1); Nrho=Nspin^2;

Mseg=zeros(n,3); BmolSeg=zeros(n,3); BexSeg=zeros(n,3); BeffSeg=zeros(n,3);
spinSeg=zeros(n,3); BmolPerpSeg=zeros(n,1); BexPerpSeg=zeros(n,1);
BeffPerpSeg=zeros(n,1); spinPerpSeg=zeros(n,1); dMmagSeg=zeros(n,1);
dMzSeg=zeros(n,1); tiltSeg=zeros(n,1); traceSeg=zeros(n,1);
hermSeg=zeros(n,1); MnormSeg=zeros(n,1);

for it=1:n
    rhoRaw=reshape(YSeg(it,1:Nrho).',Nspin,Nspin);
    rho=0.5*(rhoRaw+rhoRaw');
    Mraw=real(YSeg(it,Nrho+1:Nrho+3)); Mraw=Mraw(:).';
    Mnorm=norm(Mraw); M=Mraw/max(Mnorm,1e-15);

    J_eV=compute_bond_current_general(rho,model);
    Bmol_mT=apply_BS_kernel(model.BS_kernel_Rstar,J_eV);
    Bmol_T=1e-3*Bmol_mT;
    Bex_T=compute_Bex_fast(rho,J_site_eV,model,p);
    Beff_T=Bmol_T+Bex_T;

    obs=density_observables_from_rho(rho,p);
    sTot=[sum(obs.m_x),sum(obs.m_y),sum(obs.m_z)];

    BmolPerpVec=Bmol_T-dot(Bmol_T,M)*M;
    BexPerpVec=Bex_T-dot(Bex_T,M)*M;
    BeffPerpVec=Beff_T-dot(Beff_T,M)*M;
    spinPerpVec=sTot-dot(sTot,M)*M;

    % Full vector is the actual driving field here.
    B_LLG_T=Beff_T;
    precession=cross(M,B_LLG_T);
    damping=cross(M,cross(M,B_LLG_T));
    dM=p.gamma0_fs_T*(precession-p.lambda*damping);

    Mseg(it,:)=M; BmolSeg(it,:)=Bmol_mT; BexSeg(it,:)=Bex_T;
    BeffSeg(it,:)=Beff_T; spinSeg(it,:)=sTot;
    BmolPerpSeg(it)=norm(BmolPerpVec); BexPerpSeg(it)=norm(BexPerpVec);
    BeffPerpSeg(it)=norm(BeffPerpVec); spinPerpSeg(it)=norm(spinPerpVec);
    dMmagSeg(it)=norm(dM); dMzSeg(it)=dM(3);
    tiltSeg(it)=acosd(max(-1,min(1,M(1))));
    traceSeg(it)=abs(real(trace(rho))-1);
    hermSeg(it)=norm(rhoRaw-rhoRaw','fro');
    MnormSeg(it)=abs(Mnorm-1);
end

tObs=[tObs;tSeg]; MObs=[MObs;Mseg]; BmolObs=[BmolObs;BmolSeg];
BexObs=[BexObs;BexSeg]; BeffObs=[BeffObs;BeffSeg]; spinObs=[spinObs;spinSeg];
BmolPerpObs=[BmolPerpObs;BmolPerpSeg]; BexPerpObs=[BexPerpObs;BexPerpSeg];
BeffPerpObs=[BeffPerpObs;BeffPerpSeg]; spinPerpObs=[spinPerpObs;spinPerpSeg];
dMmagObs=[dMmagObs;dMmagSeg]; dMzObs=[dMzObs;dMzSeg]; tiltObs=[tiltObs;tiltSeg];
traceObs=[traceObs;traceSeg]; hermObs=[hermObs;hermSeg]; MnormObs=[MnormObs;MnormSeg];
lastState=YSeg(end,:).';
end

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
if isfield(E,'J0_meV')
    tagJ = safe_number_tag(E.J0_meV);
else
    tagJ = 'NA';
end

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
    sprintf('J0_%s_SOC_%s_meV_EVERY_SOLVER_STEP.csv',tagJ,tagS)));

Tf = table( ...
    (1:numel(E.finalObs.m_z)).', ...
    E.finalObs.n_up(:),E.finalObs.n_down(:), ...
    E.finalObs.m_x(:),E.finalObs.m_y(:),E.finalObs.m_z(:), ...
    'VariableNames',{'site','n_up','n_down','m_x','m_y','m_z'});

writetable(Tf,fullfile(outputFolder, ...
    sprintf('J0_%s_SOC_%s_meV_FINAL_SITE_SPIN.csv',tagJ,tagS)));

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
