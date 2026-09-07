function soc_zonly_j0()
% HELIX_J0ZERO_SOC_BOARD_NN_5MEV_BZONLY_10PS
% =========================================================================
% CLEAN SOC-ONLY CONTROL FOR THE 13-SITE HELIX MODEL
%
% QUESTION TESTED
% -------------------------------------------------------------------------
% Does molecular SOC, BY ITSELF and with NO exchange coupling, change the
% molecular Biot-Savart B_z(t) enough to reduce temporal sign cancellation
% and create a larger out-of-plane surface-moment component M_z?
%
% WHAT IS HELD FIXED
% -------------------------------------------------------------------------
%   J0                = 0 meV  -> NO exchange Hamiltonian, NO B_ex
%   M(0)              = [1 0 0]
%   lambda_LLG        = 0.6
%   tau_rel           = 500 fs
%   runtime           = 10 ps
%   same 13-site helix, contact potential, Lindblad bath, Biot-Savart,
%   ode45 tolerances, and BOARD-NN geometry-dependent SOC construction.
%
% ONLY PARAMETER SCANNED
% -------------------------------------------------------------------------
%   gamma*t = [0, +0.1, +5, -5] meV
%
% IMPORTANT LLG DIAGNOSTIC RESTRICTION
% -------------------------------------------------------------------------
% The ELECTRONIC dynamics uses the full Hamiltonian with SOC.
% The Biot-Savart field B_mol=(Bx,By,Bz) is computed normally and saved.
% However, ONLY its signed z component is passed to the LLG:
%
%       B_LLG(t) = [0, 0, B_mol,z(t)] .
%
% This is intentional so the comparison with the earlier J0=0 Bz-only
% model is apples-to-apples.  We do NOT use |Bz| in the dynamics here.
% |Bz| is used only afterwards as a cancellation/headroom diagnostic.
%
% PRIMARY OUTPUTS
% -------------------------------------------------------------------------
%   1) Bz(t) for BOARD-NN SOC = 0,+0.1,+5,-5 meV
%   2) signed integral int Bz dt
%   3) absolute integral int |Bz| dt
%   4) cancellation ratio C = |int Bz dt| / int |Bz| dt
%   5) Mz(t)
%   6) exact z-only LLG estimate Mz=tanh(lambda*gamma*int Bz dt)
%   7) electronic spin components, saved for mechanistic inspection
%
% INTERPRETATION
% -------------------------------------------------------------------------
% If BOARD-NN SOC changes Bz(t), its signed area, or the cancellation ratio at J0=0,
% then SOC already affects the molecular-current/Biot-Savart channel.
% If SOC has little effect here but a strong effect at J0=150 meV, then the
% large SOC response in the exchange model is mainly transferred through
% SOC -> transverse electronic spin -> exchange field -> LLG torque.
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

p.M0 = [1 0 0];
p.M0 = p.M0/norm(p.M0);

p.gamma0_fs_T = 1.760859e-4;
p.lambda = 0.60;

% This is the defining control of the present run.
p.J0_meV = 0;
p.lambdaSOCList_meV = [0 0.1 5 -5];

p.kB_eV_K = 8.617333262e-5;
p.T_K = 300;
p.beta_eV_inv = 1/(p.kB_eV_K*p.T_K);
p.tau_rel_fs = 500;
p.Gamma0_per_fs = 1/p.tau_rel_fs;
p.hbar_eV_fs = 0.6582119569;
p.muB_eV_T = 5.7883818060e-5;

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

p.useParallelIfAvailable = true;
p.requestedWorkers = min(3,numel(p.lambdaSOCList_meV));

%% ========================================================================
% 2. AUDIT
% =========================================================================

assert(isequal(p.M0,[1 0 0]),'M(0) must be +X.');
assert(p.J0_meV==0,'This control must have J0=0.');
assert(abs(p.lambda-0.6)<1e-14,'lambda_LLG must be 0.6.');
assert(any(p.lambdaSOCList_meV==0),'SOC=0 control must be present.');

fprintf('\n============================================================\n');
fprintf('J0=0 | BOARD-NN SOC | SIGNED Bz-ONLY LLG CONTROL\n');
fprintf('============================================================\n');
fprintf('J0               = 0 meV (exchange OFF)\n');
fprintf('SOC scan          = '); fprintf('%g ',p.lambdaSOCList_meV); fprintf('meV\n');
fprintf('M(0)              = [%g %g %g]\n',p.M0);
fprintf('lambda_LLG        = %.3f\n',p.lambda);
fprintf('tau_rel           = %.1f fs\n',p.tau_rel_fs);
fprintf('physical time     = %.2f ps\n',p.tEnd_fs/1000);
fprintf('LLG field         = [0,0,B_mol,z] ONLY (SIGNED)\n');
fprintf('electronic Zeeman = OFF\n');
fprintf('SOC model         = BOARD-NN: i*(gamma*t)*(S . Rhat_ij), gamma*t scanned in meV\n');
fprintf('solver            = ode45, Refine=1, no downsampling\n');
fprintf('============================================================\n');

%% ========================================================================
% 3. GEOMETRY + STATIC HAMILTONIAN + SOC UNIT MATRIX
% =========================================================================

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

[model.Hsoc_unit,model.soc_Rhat,model.soc_nn_pairs] = ...
    build_board_nn_soc_unit(model);

assert(norm(model.Hsoc_unit-model.Hsoc_unit','fro')<1e-12, ...
    'H_SOC unit matrix is not Hermitian.');

rhatNorm = vecnorm(model.soc_Rhat,2,2);
fprintf('\nSOC geometry audit: NN links=%d | |Rhat| %.6g ... %.6g\n', ...
    size(model.soc_nn_pairs,1),min(rhatNorm),max(rhatNorm));

model = precompute_fast_indices(model);
model.BS_kernel_Rstar = precompute_BS_kernel_at_points( ...
    Ri,bonds,p.Rstar_A,p.Nseg_BS,p.min_distance_A);

%% ========================================================================
% 4. INITIAL STATE + SAME LINDBLAD BATH FOR ALL SOC CASES
% =========================================================================

[rho0,initialDiagnostic] = build_gas_ground_initial_state(H_hop);
bath = build_fast_lindblad_bath(H_surface_orbital,p);

dr0 = apply_fast_lindblad(rho0,bath);
assert(abs(trace(dr0))<1e-11,'Lindblad trace check failed.');

fprintf('Initial rho: trace err=%.3e | Herm err=%.3e | total Mz=%.3e\n', ...
    initialDiagnostic.traceError,initialDiagnostic.hermiticityError, ...
    initialDiagnostic.initialTotalMz);

%% ========================================================================
% 5. VALIDATE PHYSICAL NN BACKBONE CURRENT AT SOC=0
% =========================================================================

case0 = build_soc_case_model(model,0);
Jold = compute_bond_current_scalar_reference(rho0,model);
Jnew = compute_bond_current_general(rho0,case0);
assert(norm(Jold-Jnew)<1e-13,'SOC-capable current failed SOC=0 baseline.');

%% ========================================================================
% 6. OUTPUT FOLDER
% =========================================================================

downloadsDir = get_downloads_directory();
runTag = datestr(now,'yyyymmdd_HHMMSS');
outputFolder = fullfile(downloadsDir, ...
    ['J0ZERO_SOC_BOARD_NN_5MEV_BZONLY_LAMBDA06_10PS_' runTag]);
if ~exist(outputFolder,'dir'), mkdir(outputFolder); end
fprintf('\nOutput folder:\n%s\n',outputFolder);

%% ========================================================================
% 7. PARALLEL AVAILABILITY
% =========================================================================

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

%% ========================================================================
% 8. RUN SOC CASES
% =========================================================================

nSOC = numel(p.lambdaSOCList_meV);
runs = cell(nSOC,1);

if useParallel
    parfor is=1:nSOC
        runs{is} = run_one_bzonly_soc_case(rho0,p.lambdaSOCList_meV(is),model,bath,p);
    end
else
    for is=1:nSOC
        runs{is} = run_one_bzonly_soc_case(rho0,p.lambdaSOCList_meV(is),model,bath,p);
    end
end

%% ========================================================================
% 9. CANCELLATION / Mz SUMMARY
% =========================================================================

soc_col = zeros(nSOC,1);
Mz_end = zeros(nSOC,1);
maxAbsMz = zeros(nSOC,1);
maxTilt_deg = zeros(nSOC,1);
maxAbsBz_mT = zeros(nSOC,1);
signedArea_Tfs = zeros(nSOC,1);
absArea_Tfs = zeros(nSOC,1);
cancelRatio_end = zeros(nSOC,1);
Mz_exact_end = zeros(nSOC,1);
Mz_abs_headroom_end = zeros(nSOC,1);
maxElectronicSpinZ = zeros(nSOC,1);
runtime_min = zeros(nSOC,1);
nSolverPoints = zeros(nSOC,1);

for is=1:nSOC
    E = runs{is};
    soc_col(is) = E.lambdaSOC_meV;
    t = E.t_fs(:);
    Bz_T = 1e-3*E.B_mol_mT(:,3);

    A_signed = cumtrapz(t,Bz_T);
    A_abs = cumtrapz(t,abs(Bz_T));
    C = abs(A_signed)./max(A_abs,eps);
    C(1) = NaN;

    E.cumBz_signed_Tfs = A_signed;
    E.cumBz_abs_Tfs = A_abs;
    E.cancellationRatio = C;
    E.Mz_exact_signed = tanh(p.lambda*p.gamma0_fs_T*A_signed);
    E.Mz_abs_headroom = tanh(p.lambda*p.gamma0_fs_T*A_abs);
    runs{is} = E;

    Mz_end(is) = E.M(end,3);
    maxAbsMz(is) = max(abs(E.M(:,3)));
    maxTilt_deg(is) = max(E.tilt_deg);
    maxAbsBz_mT(is) = max(abs(E.B_mol_mT(:,3)));
    signedArea_Tfs(is) = A_signed(end);
    absArea_Tfs(is) = A_abs(end);
    cancelRatio_end(is) = abs(A_signed(end))/max(A_abs(end),eps);
    Mz_exact_end(is) = E.Mz_exact_signed(end);
    Mz_abs_headroom_end(is) = E.Mz_abs_headroom(end);
    maxElectronicSpinZ(is) = max(abs(E.totalElectronicSpin(:,3)));
    runtime_min(is) = E.runtime_min;
    nSolverPoints(is) = numel(E.t_fs);

    save(fullfile(outputFolder,sprintf('SOC_%s_meV_COMPACT.mat', ...
        safe_number_tag(E.lambdaSOC_meV))),'E','-v7.3');
    write_bz_case_csv(E,outputFolder);
end

summaryTable = table( ...
    soc_col,Mz_end,maxAbsMz,maxTilt_deg,maxAbsBz_mT, ...
    signedArea_Tfs,absArea_Tfs,cancelRatio_end, ...
    Mz_exact_end,Mz_abs_headroom_end,maxElectronicSpinZ, ...
    nSolverPoints,runtime_min, ...
    'VariableNames',{ ...
    'lambda_SOC_meV','Mz_end','max_abs_Mz','max_tilt_deg', ...
    'max_abs_Bz_mT','signed_Bz_area_Tfs','abs_Bz_area_Tfs', ...
    'cancellation_ratio_end','Mz_exact_signed_end', ...
    'Mz_abs_headroom_end','max_abs_electronic_spin_z', ...
    'n_solver_points','runtime_min'});

writetable(summaryTable,fullfile(outputFolder,'J0ZERO_BOARD_NN_SOC_BZONLY_SUMMARY.csv'));
save(fullfile(outputFolder,'COMPLETE_J0ZERO_BOARD_NN_SOC_BZONLY.mat'), ...
    'runs','summaryTable','p','model','initialDiagnostic','-v7.3');

fprintf('\n============================================================\n');
fprintf('J0=0 BOARD-NN SOC Bz-ONLY SUMMARY\n');
fprintf('============================================================\n');
disp(summaryTable);

for is=1:nSOC
    fprintf(['SOC=%+g meV | Mz(end)=%+.6e | intBz=%+.6e Tfs | ' ...
             'int|Bz|=%.6e Tfs | C_end=%.4g\n'], ...
        soc_col(is),Mz_end(is),signedArea_Tfs(is),absArea_Tfs(is), ...
        cancelRatio_end(is));
end

%% ========================================================================
% 10. FIGURES
% =========================================================================

labels = arrayfun(@(x)sprintf('SOC = %+g meV',x),soc_col,'UniformOutput',false);

% FIG 0: geometry / BOARD-NN bond directions
fig0 = figure('Color','w','Position',[60 60 1200 850]);
plot3(model.Ri(:,1),model.Ri(:,2),model.Ri(:,3),'-o','LineWidth',1.2); hold on;
for ii=1:size(model.soc_nn_pairs,1)
    i = model.soc_nn_pairs(ii,1); j = model.soc_nn_pairs(ii,2);
    rmid = 0.5*(model.Ri(i,:)+model.Ri(j,:)); rh = model.soc_Rhat(ii,:);
    quiver3(rmid(1),rmid(2),rmid(3),rh(1),rh(2),rh(3),0.55,'LineWidth',1.0);
end
grid on; axis equal; xlabel('x [A]'); ylabel('y [A]'); zlabel('z [A]');
title('Helix geometry and BOARD-NN SOC bond directions Rhat_{ij}'); view(35,25);
save_figure(fig0,outputFolder,'FIG00_BOARD_NN_SOC_GEOMETRY.png');

% FIG 1: Bz(t), early + full
fig1 = figure('Color','w','Position',[80 80 1450 850]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile; hold on;
for is=1:nSOC
    E=runs{is}; plot(E.t_fs/1000,E.B_mol_mT(:,3),'LineWidth',1.2,'DisplayName',labels{is});
end
grid on; xlim([0 1]); ylabel('B_z [mT]');
title('Early-time molecular B_z: does SOC alter the sign oscillations?'); legend('Location','best');
nexttile; hold on;
for is=1:nSOC
    E=runs{is}; plot(E.t_fs/1000,E.B_mol_mT(:,3),'LineWidth',1.2,'DisplayName',labels{is});
end
grid on; xlim([0 p.tEnd_fs/1000]); xlabel('time [ps]'); ylabel('B_z [mT]');
title('Full 10 ps molecular B_z'); legend('Location','best');
save_figure(fig1,outputFolder,'FIG01_Bz_SOC_SCAN.png');

% FIG 2: signed and absolute accumulated area for each SOC case
fig2 = figure('Color','w','Position',[100 70 1500 930]);
tiledlayout(nSOC,1,'TileSpacing','compact','Padding','compact');
for is=1:nSOC
    E=runs{is}; nexttile;
    plot(E.t_fs/1000,E.cumBz_signed_Tfs,'LineWidth',1.25,'DisplayName','int B_z dt'); hold on;
    plot(E.t_fs/1000,E.cumBz_abs_Tfs,'LineWidth',1.25,'DisplayName','int |B_z| dt');
    grid on; ylabel('T fs'); title(labels{is}); legend('Location','best');
    if is==nSOC, xlabel('time [ps]'); end
end
sgtitle('SOC effect on signed versus absolute accumulated B_z');
save_figure(fig2,outputFolder,'FIG02_ACCUMULATED_Bz_SOC_SCAN.png');

% FIG 3: cancellation ratio C(t)
fig3 = figure('Color','w','Position',[120 100 1400 760]); hold on;
for is=1:nSOC
    E=runs{is}; plot(E.t_fs/1000,E.cancellationRatio,'LineWidth',1.25,'DisplayName',labels{is});
end
grid on; xlabel('time [ps]'); ylabel('|int B_z dt| / int |B_z| dt');
title('Temporal cancellation ratio: does SOC rectify B_z?'); legend('Location','best');
xlim([0 p.tEnd_fs/1000]);
save_figure(fig3,outputFolder,'FIG03_CANCELLATION_RATIO_SOC_SCAN.png');

% FIG 4: Mz(t)
fig4 = figure('Color','w','Position',[140 100 1400 760]); hold on;
for is=1:nSOC
    E=runs{is}; plot(E.t_fs/1000,E.M(:,3),'LineWidth',1.35,'DisplayName',labels{is});
end
grid on; xlabel('time [ps]'); ylabel('M_z');
title('Out-of-plane surface-moment component | J0=0 | signed B_z-only LLG');
legend('Location','best'); xlim([0 p.tEnd_fs/1000]);
save_figure(fig4,outputFolder,'FIG04_Mz_SOC_SCAN.png');

% FIG 5: actual LLG Mz vs exact z-only integral formula
fig5 = figure('Color','w','Position',[160 60 1500 930]);
tiledlayout(nSOC,1,'TileSpacing','compact','Padding','compact');
for is=1:nSOC
    E=runs{is}; nexttile;
    plot(E.t_fs/1000,E.M(:,3),'LineWidth',1.25,'DisplayName','LLG M_z'); hold on;
    plot(E.t_fs/1000,E.Mz_exact_signed,'--','LineWidth',1.25, ...
        'DisplayName','tanh(lambda gamma int B_z dt)');
    grid on; ylabel('M_z'); title(labels{is}); legend('Location','best');
    if is==nSOC, xlabel('time [ps]'); end
end
sgtitle('Direct consistency test: z-only LLG versus accumulated signed B_z');
save_figure(fig5,outputFolder,'FIG05_LLG_VS_INTEGRAL_SOC_SCAN.png');

% FIG 6: final metrics versus SOC
[socPlot,ord] = sort(soc_col);
fig6 = figure('Color','w','Position',[180 80 1450 900]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile; plot(socPlot,Mz_end(ord),'-o','LineWidth',1.4); yline(0,'--'); grid on;
xlabel('lambda_SOC [meV]'); ylabel('M_z(end)'); title('Final out-of-plane component');
nexttile; plot(socPlot,signedArea_Tfs(ord),'-o','LineWidth',1.4); yline(0,'--'); grid on;
xlabel('lambda_SOC [meV]'); ylabel('int B_z dt [T fs]'); title('Signed accumulated molecular B_z');
nexttile; plot(socPlot,cancelRatio_end(ord),'-o','LineWidth',1.4); grid on;
xlabel('lambda_SOC [meV]'); ylabel('C_{end}'); title('Remaining signed fraction');
nexttile; plot(socPlot,maxAbsBz_mT(ord),'-o','LineWidth',1.4); grid on;
xlabel('lambda_SOC [meV]'); ylabel('max |B_z| [mT]'); title('Instantaneous B_z amplitude');
sgtitle('What does SOC change at J0=0?');
save_figure(fig6,outputFolder,'FIG06_SOC_EFFECT_SUMMARY_J0ZERO.png');

% FIG 7: actual signed Mz versus artificial |Bz| headroom estimate
fig7 = figure('Color','w','Position',[200 100 1400 760]); hold on;
for is=1:nSOC
    E=runs{is};
    plot(E.t_fs/1000,E.M(:,3),'LineWidth',1.2,'DisplayName',[labels{is} ' actual signed']);
    plot(E.t_fs/1000,E.Mz_abs_headroom,'--','LineWidth',1.1, ...
        'DisplayName',[labels{is} ' |B_z| headroom']);
end
grid on; xlabel('time [ps]'); ylabel('M_z');
title('|B_z| shown ONLY as cancellation headroom diagnostic; not used in dynamics');
legend('Location','best'); xlim([0 p.tEnd_fs/1000]);
save_figure(fig7,outputFolder,'FIG07_SIGNED_VS_ABS_HEADROOM.png');

%% ========================================================================
% 11. MODEL NOTE
% =========================================================================

fid=fopen(fullfile(outputFolder,'MODEL_NOTE.txt'),'w');
if fid>=0
    fprintf(fid,'J0=0 BOARD-NN SOC Bz-only control\n');
    fprintf(fid,'SOC scan: 0,+0.1,+5,-5 meV.\n');
    fprintf(fid,'No exchange and no electronic Zeeman.\n');
    fprintf(fid,'Electronic dynamics uses BOARD-NN SOC: Hsoc_ij = i*(gamma*t)*(S.Rhat_ij), with S=sigma/2.\n');
    fprintf(fid,'Here gamma*t is the scanned energy parameter; +5 meV follows the board value.\n');
    fprintf(fid,'LLG receives only signed [0,0,B_mol,z].\n');
    fprintf(fid,'|Bz| is post-processing only; it is NOT used in the physical run.\n');
    fclose(fid);
end

fprintf('\n============================================================\n');
fprintf('RUN COMPLETE\nOutput folder:\n%s\n',outputFolder);
fprintf('Total wall time = %.2f min\n',toc(totalTimer)/60);
fprintf('============================================================\n');

end


%% =========================================================================
% ONE SOC CASE: J0=0, SIGNED Bz-ONLY LLG
% =========================================================================

function E = run_one_bzonly_soc_case(rho0,soc_meV,model,bath,p)

caseTimer=tic;
fprintf('\n------------------------------------------------------------\n');
fprintf('START J0=0 | SOC=%+g meV | Bz-only LLG\n',soc_meV);
fprintf('------------------------------------------------------------\n');

caseModel=build_soc_case_model(model,soc_meV);
Nspin=size(caseModel.H_static_spin,1);
Nrho=Nspin^2;
rhs=@(~,y) bzonly_soc_rhs(y,caseModel,bath,p);
y0=[rho0(:);p.M0(:)];

% Every accepted ode45 point is retained.
tObs=[]; MObs=[]; BmolObs=[]; BllgObs=[]; spinObs=[];
dMmagObs=[]; dMzObs=[]; tiltObs=[]; traceObs=[]; hermObs=[]; MnormObs=[];

optsEarly=odeset('InitialStep',p.InitialStep_fs, ...
    'MaxStep',p.MaxStepEarly_fs,'RelTol',p.RelTol,'AbsTol',p.AbsTol, ...
    'Refine',1,'Stats','off');
[tSeg,YSeg]=ode45(rhs,[p.tStart_fs p.earlyEnd_fs],y0,optsEarly);

[tObs,MObs,BmolObs,BllgObs,spinObs,dMmagObs,dMzObs,tiltObs, ...
 traceObs,hermObs,MnormObs,lastState] = append_bzonly_segment( ...
 tSeg,YSeg,false,tObs,MObs,BmolObs,BllgObs,spinObs,dMmagObs,dMzObs, ...
 tiltObs,traceObs,hermObs,MnormObs,caseModel,p);

fprintf('  SOC=%+g: reached %.3f ps | points=%d\n', ...
    soc_meV,p.earlyEnd_fs/1000,numel(tSeg));

optsLate=odeset('RelTol',p.RelTol,'AbsTol',p.AbsTol,'Refine',1,'Stats','off');
currentTime=p.earlyEnd_fs;
while currentTime<p.tEnd_fs-1e-12
    nextEnd=min(currentTime+p.chunk_fs,p.tEnd_fs);
    segTimer=tic;
    [tSeg,YSeg]=ode45(rhs,[currentTime nextEnd],lastState,optsLate);
    [tObs,MObs,BmolObs,BllgObs,spinObs,dMmagObs,dMzObs,tiltObs, ...
     traceObs,hermObs,MnormObs,lastState] = append_bzonly_segment( ...
     tSeg,YSeg,true,tObs,MObs,BmolObs,BllgObs,spinObs,dMmagObs,dMzObs, ...
     tiltObs,traceObs,hermObs,MnormObs,caseModel,p);
    currentTime=nextEnd;
    dt=diff(tSeg); if isempty(dt), dtMed=NaN; else, dtMed=median(dt); end
    fprintf('  SOC=%+g: reached %5.2f ps | points=%7d | median dt=%g fs | %.2f min\n', ...
        soc_meV,currentTime/1000,numel(tSeg),dtMed,toc(segTimer)/60);
end

rhoFinal=reshape(lastState(1:Nrho),Nspin,Nspin);
rhoFinal=0.5*(rhoFinal+rhoFinal');

E.lambdaSOC_meV=soc_meV;
E.t_fs=tObs;
E.M=MObs;
E.B_mol_mT=BmolObs;
E.B_LLG_T=BllgObs;
E.totalElectronicSpin=spinObs;
E.dM_mag_per_fs=dMmagObs;
E.dMz_per_fs=dMzObs;
E.tilt_deg=tiltObs;
E.traceError=traceObs;
E.hermiticityError=hermObs;
E.MnormError=MnormObs;
E.finalRho=rhoFinal;
E.finalObs=density_observables_from_rho(rhoFinal,p);
E.runtime_min=toc(caseTimer)/60;

fprintf('COMPLETE SOC=%+g | runtime %.2f min | final Mz=%+.6e | max|Bz|=%.6f mT\n', ...
    soc_meV,E.runtime_min,E.M(end,3),max(abs(E.B_mol_mT(:,3))));

end


%% =========================================================================
% RHS: ELECTRONS WITH SOC; SURFACE MOMENT SEES ONLY SIGNED Bz
% =========================================================================

function dy = bzonly_soc_rhs(y,model,bath,p)

Nspin=size(model.H_static_spin,1);
Nrho=Nspin^2;
rho=reshape(y(1:Nrho),Nspin,Nspin);
rho=0.5*(rho+rho');

M=real(y(Nrho+1:Nrho+3));
M=M(:).'; M=M/max(norm(M),1e-15);

% Physical NN currents. SOC changes rho, hence it can alter these currents.
J_eV=compute_bond_current_general(rho,model);
B_mol_mT=apply_BS_kernel(model.BS_kernel_Rstar,J_eV);
B_mol_T=1e-3*B_mol_mT;

% J0=0: NO H_MS and NO exchange feedback.
H_total=model.H_static_spin;
dr_H=-1i/p.hbar_eV_fs*(H_total*rho-rho*H_total);
dr_D=apply_fast_lindblad(rho,bath);
drho=dr_H+dr_D;

% IMPORTANT: only the signed z component drives LLG.
B_LLG_T=[0 0 B_mol_T(3)];
precession=cross(M,B_LLG_T);
damping=cross(M,cross(M,B_LLG_T));
dM=p.gamma0_fs_T*(precession-p.lambda*damping);

dy=[drho(:);dM(:)];

end


%% =========================================================================
% APPEND ALL ACCEPTED SOLVER POINTS
% =========================================================================

function [tObs,MObs,BmolObs,BllgObs,spinObs,dMmagObs,dMzObs,tiltObs, ...
          traceObs,hermObs,MnormObs,lastState] = append_bzonly_segment( ...
    tSeg,YSeg,dropFirst,tObs,MObs,BmolObs,BllgObs,spinObs,dMmagObs,dMzObs, ...
    tiltObs,traceObs,hermObs,MnormObs,model,p)

if dropFirst
    tSeg=tSeg(2:end); YSeg=YSeg(2:end,:);
end

n=numel(tSeg);
Nspin=size(model.H_static_spin,1);
Nrho=Nspin^2;

Mseg=zeros(n,3); BmolSeg=zeros(n,3); BllgSeg=zeros(n,3); spinSeg=zeros(n,3);
dMmagSeg=zeros(n,1); dMzSeg=zeros(n,1); tiltSeg=zeros(n,1);
traceSeg=zeros(n,1); hermSeg=zeros(n,1); MnormSeg=zeros(n,1);

for it=1:n
    rhoRaw=reshape(YSeg(it,1:Nrho).',Nspin,Nspin);
    rho=0.5*(rhoRaw+rhoRaw');
    Mraw=real(YSeg(it,Nrho+1:Nrho+3)); Mraw=Mraw(:).';
    Mnorm=norm(Mraw); M=Mraw/max(Mnorm,1e-15);

    J_eV=compute_bond_current_general(rho,model);
    Bmol_mT=apply_BS_kernel(model.BS_kernel_Rstar,J_eV);
    Bmol_T=1e-3*Bmol_mT;
    Bllg=[0 0 Bmol_T(3)];

    obs=density_observables_from_rho(rho,p);
    sTot=[sum(obs.m_x),sum(obs.m_y),sum(obs.m_z)];

    precession=cross(M,Bllg);
    damping=cross(M,cross(M,Bllg));
    dM=p.gamma0_fs_T*(precession-p.lambda*damping);

    Mseg(it,:)=M;
    BmolSeg(it,:)=Bmol_mT;
    BllgSeg(it,:)=Bllg;
    spinSeg(it,:)=sTot;
    dMmagSeg(it)=norm(dM);
    dMzSeg(it)=dM(3);
    tiltSeg(it)=acosd(max(-1,min(1,M(1))));
    traceSeg(it)=abs(real(trace(rho))-1);
    hermSeg(it)=norm(rhoRaw-rhoRaw','fro');
    MnormSeg(it)=abs(Mnorm-1);
end

tObs=[tObs;tSeg]; %#ok<AGROW>
MObs=[MObs;Mseg]; %#ok<AGROW>
BmolObs=[BmolObs;BmolSeg]; %#ok<AGROW>
BllgObs=[BllgObs;BllgSeg]; %#ok<AGROW>
spinObs=[spinObs;spinSeg]; %#ok<AGROW>
dMmagObs=[dMmagObs;dMmagSeg]; %#ok<AGROW>
dMzObs=[dMzObs;dMzSeg]; %#ok<AGROW>
tiltObs=[tiltObs;tiltSeg]; %#ok<AGROW>
traceObs=[traceObs;traceSeg]; %#ok<AGROW>
hermObs=[hermObs;hermSeg]; %#ok<AGROW>
MnormObs=[MnormObs;MnormSeg]; %#ok<AGROW>
lastState=YSeg(end,:).';

end


%% =========================================================================
% CSV WRITER FOR THIS CONTROL
% =========================================================================

function write_bz_case_csv(E,outputFolder)

tagS=safe_number_tag(E.lambdaSOC_meV);
T=table( ...
    E.t_fs(:),E.t_fs(:)/1000, ...
    E.M(:,1),E.M(:,2),E.M(:,3), ...
    E.B_mol_mT(:,1),E.B_mol_mT(:,2),E.B_mol_mT(:,3), ...
    E.B_LLG_T(:,1),E.B_LLG_T(:,2),E.B_LLG_T(:,3), ...
    E.totalElectronicSpin(:,1),E.totalElectronicSpin(:,2),E.totalElectronicSpin(:,3), ...
    E.cumBz_signed_Tfs(:),E.cumBz_abs_Tfs(:),E.cancellationRatio(:), ...
    E.Mz_exact_signed(:),E.Mz_abs_headroom(:), ...
    E.dM_mag_per_fs(:),E.dMz_per_fs(:),E.tilt_deg(:), ...
    E.traceError(:),E.hermiticityError(:),E.MnormError(:), ...
    'VariableNames',{ ...
    'time_fs','time_ps','Mx','My','Mz', ...
    'Bmol_x_mT','Bmol_y_mT','Bmol_z_mT', ...
    'BLLG_x_T','BLLG_y_T','BLLG_z_T', ...
    'electronic_mx_total','electronic_my_total','electronic_mz_total', ...
    'cum_Bz_signed_Tfs','cum_Bz_abs_Tfs','cancellation_ratio', ...
    'Mz_exact_signed','Mz_abs_headroom', ...
    'dM_mag_per_fs','dMz_per_fs','tilt_deg', ...
    'trace_error','hermiticity_error','Mnorm_error'});

writetable(T,fullfile(outputFolder, ...
    sprintf('SOC_%s_meV_EVERY_SOLVER_STEP.csv',tagS)));
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

% Precompute the FULL physical NN hopping blocks, INCLUDING BOARD-NN SOC.
% Because SOC lives on the same NN bonds, the charge-current operator on each
% physical bond must use the complete 2x2 bond block H_ij.
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
% BOARD-NN SOC is part of the physical NN hopping block H_ij, so the
% current uses the full 2x2 spin-dependent bond Hamiltonian. For SOC=0,
% H_ij=-t I and this reduces exactly to
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
