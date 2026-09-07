function aniso_threshold()
% =========================================================================
% STAGE 3: Ku THRESHOLD + MECHANISM TEST | FULL XYZ LLG
%
% Goal:
%   1) Resolve the Ku threshold between 0.25 and 0.50 meV.
%   2) Add SOC=0 as a control to test whether Ku can drive reorientation
%      without BOARD-NN SOC.
%   3) Track Mx, My, Mz and Bmol/Bex/Baniso/Beff near the transition.
%   4) Quantify flip time and temporal cancellation.
%
% Fixed:
%   J0 = 150 meV, M(0)=+X, lambda_LLG=0.60, tau_rel=500 fs, 10 ps,
%   same helix, bath, Biot-Savart, exchange feedback and two-stage ode45.
%
% Scanned:
%   SOC = 0, +5, -5 meV
%   Ku  = 0.25, 0.30, 0.35, 0.40, 0.45, 0.50 meV
%
% Anisotropy field added to the LLG:
%   B_aniso(M) = (2*Ku/mu_s) * (M.e_easy) * e_easy
% with e_easy = +Z. Stage 3 focuses only on the transition window.
%
% Key diagnostic:
%   C_torque = | integral dMz/dt dt | / integral |dMz/dt| dt
%   C_torque -> 1 : little temporal cancellation
%   C_torque -> 0 : strong temporal cancellation
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

p.J0_meV = 150;
p.socList_meV = [0 5 -5];
p.r0_A = 3.0;

% Ku values are stored in eV. Stage 3 resolves the 0.25-0.50 meV threshold window.
p.KuList_eV = 1e-3*[0.25 0.30 0.35 0.40 0.45 0.50];
p.easyAxis  = [0 0 1];

p.kB_eV_K = 8.617333262e-5;
p.T_K = 300;
p.beta_eV_inv = 1/(p.kB_eV_K*p.T_K);
p.tau_rel_fs = 500;
p.Gamma0_per_fs = 1/p.tau_rel_fs;
p.hbar_eV_fs = 0.6582119569;
p.muB_eV_T = 5.7883818060e-5;
p.surface_moment_muB = 1.0;
p.surface_moment_eV_T = p.surface_moment_muB*p.muB_eV_T;

% Same validated two-stage integration as the previous controls.
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

assert(isequal(p.M0,[1 0 0]),'M(0) must be +X.');
assert(abs(p.lambda-0.6)<1e-14,'lambda_LLG must remain 0.6.');
assert(abs(p.tau_rel_fs-500)<1e-12,'tau_rel must remain 500 fs.');
assert(all(p.KuList_eV>=0.25e-3 & p.KuList_eV<=0.50e-3),'Stage 3 Ku scan must stay in 0.25-0.50 meV.');
assert(all(ismember([0 5 -5],p.socList_meV)),'SOC list must contain 0, +5 and -5 meV.');

fprintf('\n============================================================\n');
fprintf('STAGE 3 | Ku THRESHOLD + MECHANISM | FULL XYZ LLG\n');
fprintf('============================================================\n');
fprintf('J0                 = %g meV\n',p.J0_meV);
fprintf('SOC scan           = '); fprintf('%+g ',p.socList_meV); fprintf('meV\n');
fprintf('Ku scan            = '); fprintf('%g ',1e3*p.KuList_eV); fprintf('meV\n');
fprintf('easy axis          = [%g %g %g]\n',p.easyAxis);
fprintf('M(0)               = [%g %g %g]\n',p.M0);
fprintf('lambda_LLG         = %.3f\n',p.lambda);
fprintf('tau_rel            = %.1f fs\n',p.tau_rel_fs);
fprintf('physical time      = %.2f ps\n',p.tEnd_fs/1000);
fprintf('B_aniso(M)         = (2*Ku/mu_s)*(M.e)*e\n');
fprintf('LLG field          = FULL XYZ: Bmol + Bex + Baniso\n');
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

%% 4. VALIDATE SOC-CAPABLE NN CURRENT AT SOC=0 (unchanged check) =========
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
    ['ANISOTROPY_STAGE3_THRESHOLD_MECHANISM_J150_FULLXYZ_10PS_' runTag]);
if ~exist(outputFolder,'dir'), mkdir(outputFolder); end
fprintf('\nOutput folder:\n%s\n',outputFolder);

%% 6. EXCHANGE DATA =======================================================
J_site_eV = exchange_profile(p.J0_meV,model,p);
exchange = precompute_exchange_matrices(J_site_eV);

%% 7. RUN SOC x Ku MATRIX =================================================
nSOC = numel(p.socList_meV);
nKu  = numel(p.KuList_eV);
runs = cell(nSOC,nKu);

for isoc = 1:nSOC
    soc_meV = p.socList_meV(isoc);
    fprintf('\n============================================================\n');
    fprintf('START SOC = %+g meV\n',soc_meV);
    fprintf('============================================================\n');

    % SOC is part of the electronic Hamiltonian for this whole Ku row.
    caseModelSOC = build_soc_case_model(model,soc_meV);

    for ik = 1:nKu
        Ku_eV = p.KuList_eV(ik);
        aniso = build_aniso_case_model(Ku_eV,p.easyAxis,p.surface_moment_eV_T);

        runs{isoc,ik} = run_one_case_aniso( ...
            rho0,p.J0_meV,soc_meV,aniso, ...
            J_site_eV,exchange,caseModelSOC,bath,p,'FULLXYZ');
    end
end

%% 8. BASIC RUN AUDIT ================================================
fprintf('\nStage-3 run audit: final M for the lowest Ku = %.2f meV\n',1e3*p.KuList_eV(1));
for isoc = 1:nSOC
    E0 = runs{isoc,1};
    fprintf('  SOC=%+g meV | final M=[%.6g %.6g %.6g] | max|Mz|=%.6e\n', ...
        E0.lambdaSOC_meV,E0.M(end,1),E0.M(end,2),E0.M(end,3),max(abs(E0.M(:,3))));
end

%% 9. SUMMARY TABLE + CANCELLATION METRIC =================================
nRows = nSOC*nKu;
SOC_col = zeros(nRows,1);
Ku_col_meV = zeros(nRows,1);
Mz_end = zeros(nRows,1);
max_abs_Mz = zeros(nRows,1);
torque_cancellation_ratio = zeros(nRows,1);
max_Beff_perp_mT = zeros(nRows,1);
max_Baniso_mT = zeros(nRows,1);
max_tilt_deg = zeros(nRows,1);
flip_time_ps = nan(nRows,1);
Bmol_perp_at_flip_mT = nan(nRows,1);
Bex_perp_at_flip_mT = nan(nRows,1);
Baniso_at_flip_mT = nan(nRows,1);
Beff_perp_at_flip_mT = nan(nRows,1);

row = 0;
for isoc = 1:nSOC
    for ik = 1:nKu
        row = row + 1;
        E = runs{isoc,ik};

        SOC_col(row) = E.lambdaSOC_meV;
        Ku_col_meV(row) = 1e3*E.Ku_eV;
        Mz_end(row) = E.M(end,3);
        max_abs_Mz(row) = max(abs(E.M(:,3)));

        signedTorqueArea = trapz(E.t_fs,E.dMz_per_fs);
        absTorqueArea = trapz(E.t_fs,abs(E.dMz_per_fs));
        torque_cancellation_ratio(row) = ...
            abs(signedTorqueArea)/max(absTorqueArea,eps);

        max_Beff_perp_mT(row) = 1e3*max(E.Beff_perp_T);
        max_Baniso_mT(row) = 1e3*max(vecnorm(E.B_aniso_T,2,2));
        max_tilt_deg(row) = max(E.tilt_deg);

        % Define a near-complete reorientation as |Mz| >= 0.90.
        iFlip = find(abs(E.M(:,3)) >= 0.90,1,'first');
        if ~isempty(iFlip)
            flip_time_ps(row) = E.t_fs(iFlip)/1000;
            Bmol_perp_at_flip_mT(row) = 1e3*E.Bmol_perp_T(iFlip);
            Bex_perp_at_flip_mT(row) = 1e3*E.Bex_perp_T(iFlip);
            Baniso_at_flip_mT(row) = 1e3*norm(E.B_aniso_T(iFlip,:));
            Beff_perp_at_flip_mT(row) = 1e3*E.Beff_perp_T(iFlip);
        end
    end
end

summaryTable = table( ...
    SOC_col,Ku_col_meV,Mz_end,max_abs_Mz,torque_cancellation_ratio, ...
    max_Beff_perp_mT,max_Baniso_mT,max_tilt_deg,flip_time_ps, ...
    Bmol_perp_at_flip_mT,Bex_perp_at_flip_mT,Baniso_at_flip_mT,Beff_perp_at_flip_mT, ...
    'VariableNames',{ ...
    'SOC_meV','Ku_meV','Mz_end','max_abs_Mz','torque_cancellation_ratio', ...
    'max_Beff_perp_mT','max_Baniso_mT','max_tilt_deg','flip_time_ps', ...
    'Bmol_perp_at_flip_mT','Bex_perp_at_flip_mT','Baniso_at_flip_mT','Beff_perp_at_flip_mT'});

disp(summaryTable);
writetable(summaryTable,fullfile(outputFolder,'SOC_KU_STAGE3_THRESHOLD_SUMMARY.csv'));
save(fullfile(outputFolder,'SOC_KU_STAGE3_THRESHOLD_COMPLETE.mat'), ...
    'runs','summaryTable','p','model','initialDiagnostic','-v7.3');

%% 10. FIGURES ============================================================
% A. Fine Mz(t) Ku scan, one figure for each SOC value.
for isoc = 1:nSOC
    fig = figure('Color','w','Position',[180 80 1350 800]); hold on;
    for ik = 1:nKu
        E = runs{isoc,ik};
        plot(E.t_fs/1000,E.M(:,3),'LineWidth',1.35, ...
            'DisplayName',sprintf('K_u = %.2f meV',1e3*E.Ku_eV));
    end
    yline(0,'--'); grid on;
    xlabel('time [ps]'); ylabel('M_z');
    title(sprintf('Stage 3 | SOC = %+g meV | J_0 = 150 meV | Full XYZ | fine K_u scan',p.socList_meV(isoc)));
    legend('Location','best'); xlim([0 10]);
    save_figure(fig,outputFolder,sprintf('FIG_STAGE3_MZ_SOC_%+gmeV_FINE_KU.png',p.socList_meV(isoc)));
end

% B. Final Mz versus Ku for SOC=0,+5,-5.
fig2 = figure('Color','w','Position',[220 100 1050 700]); hold on;
for isoc = 1:nSOC
    idxRow = summaryTable.SOC_meV == p.socList_meV(isoc);
    plot(summaryTable.Ku_meV(idxRow),summaryTable.Mz_end(idxRow), ...
        '-o','LineWidth',1.5,'DisplayName',sprintf('SOC = %+g meV',p.socList_meV(isoc)));
end
yline(0,'--'); grid on;
xlabel('K_u [meV]'); ylabel('M_z(10 ps)');
title('Stage 3: threshold and SOC control of the final M_z');
legend('Location','best');
save_figure(fig2,outputFolder,'FIG_STAGE3_FINAL_MZ_VS_KU.png');

% C. Temporal survival ratio. 1 = almost no cancellation.
fig3 = figure('Color','w','Position',[220 100 1050 700]); hold on;
for isoc = 1:nSOC
    idxRow = summaryTable.SOC_meV == p.socList_meV(isoc);
    plot(summaryTable.Ku_meV(idxRow),summaryTable.torque_cancellation_ratio(idxRow), ...
        '-o','LineWidth',1.5,'DisplayName',sprintf('SOC = %+g meV',p.socList_meV(isoc)));
end
grid on;
xlabel('K_u [meV]');
ylabel('|int dM_z/dt dt| / int |dM_z/dt| dt');
title('Stage 3: survival of signed M_z rotation');
legend('Location','best'); ylim([0 1]);
save_figure(fig3,outputFolder,'FIG_STAGE3_CANCELLATION_SURVIVAL_VS_KU.png');

% D. Flip time versus Ku, where flip means first |Mz| >= 0.90.
fig4 = figure('Color','w','Position',[220 100 1050 700]); hold on;
for isoc = 1:nSOC
    idxRow = summaryTable.SOC_meV == p.socList_meV(isoc);
    plot(summaryTable.Ku_meV(idxRow),summaryTable.flip_time_ps(idxRow), ...
        '-o','LineWidth',1.5,'DisplayName',sprintf('SOC = %+g meV',p.socList_meV(isoc)));
end
grid on;
xlabel('K_u [meV]'); ylabel('first time |M_z| >= 0.90 [ps]');
title('Stage 3: reorientation time versus anisotropy');
legend('Location','best');
save_figure(fig4,outputFolder,'FIG_STAGE3_FLIP_TIME_VS_KU.png');

% E. Mechanism figures at selected Ku values around the threshold.
KuDiag_meV = [0.25 0.40 0.50];
for isoc = 1:nSOC
    for q = 1:numel(KuDiag_meV)
        [~,ik] = min(abs(1e3*p.KuList_eV-KuDiag_meV(q)));
        E = runs{isoc,ik};

        % Full moment components.
        figM = figure('Color','w','Position',[180 80 1250 760]); hold on;
        plot(E.t_fs/1000,E.M(:,1),'LineWidth',1.4,'DisplayName','M_x');
        plot(E.t_fs/1000,E.M(:,2),'LineWidth',1.4,'DisplayName','M_y');
        plot(E.t_fs/1000,E.M(:,3),'LineWidth',1.6,'DisplayName','M_z');
        yline(0,'--'); grid on; xlim([0 10]); ylim([-1.05 1.05]);
        xlabel('time [ps]'); ylabel('M component');
        title(sprintf('Moment path | SOC=%+g meV | K_u=%.2f meV',E.lambdaSOC_meV,1e3*E.Ku_eV));
        legend('Location','best');
        save_figure(figM,outputFolder,sprintf('FIG_STAGE3_MCOMP_SOC_%+g_KU_%.2fmeV.png', ...
            E.lambdaSOC_meV,1e3*E.Ku_eV));

        % Competing field scales. Baniso is shown as its magnitude because it lies on the easy axis.
        figB = figure('Color','w','Position',[180 80 1250 760]); hold on;
        plot(E.t_fs/1000,1e3*E.Bmol_perp_T,'LineWidth',1.3,'DisplayName','|B_{mol,\perp}|');
        plot(E.t_fs/1000,1e3*E.Bex_perp_T,'LineWidth',1.3,'DisplayName','|B_{ex,\perp}|');
        plot(E.t_fs/1000,1e3*vecnorm(E.B_aniso_T,2,2),'LineWidth',1.5,'DisplayName','|B_{aniso}|');
        plot(E.t_fs/1000,1e3*E.Beff_perp_T,'LineWidth',1.5,'DisplayName','|B_{eff,\perp}|');
        grid on; xlim([0 10]);
        xlabel('time [ps]'); ylabel('field [mT]');
        title(sprintf('Field competition | SOC=%+g meV | K_u=%.2f meV',E.lambdaSOC_meV,1e3*E.Ku_eV));
        legend('Location','best');
        save_figure(figB,outputFolder,sprintf('FIG_STAGE3_FIELDS_SOC_%+g_KU_%.2fmeV.png', ...
            E.lambdaSOC_meV,1e3*E.Ku_eV));
    end
end

% F. At-flip field summary. Only cases that reach |Mz|>=0.90 appear.
fig5 = figure('Color','w','Position',[220 100 1100 720]); hold on;
for isoc = 1:nSOC
    idxRow = summaryTable.SOC_meV == p.socList_meV(isoc);
    plot(summaryTable.Ku_meV(idxRow),summaryTable.Baniso_at_flip_mT(idxRow), ...
        '-o','LineWidth',1.5,'DisplayName',sprintf('|B_{aniso}| at flip, SOC=%+g',p.socList_meV(isoc)));
end
grid on;
xlabel('K_u [meV]'); ylabel('|B_{aniso}| at first |M_z|>=0.90 [mT]');
title('Stage 3: anisotropy field scale at the reorientation threshold');
legend('Location','best');
save_figure(fig5,outputFolder,'FIG_STAGE3_BANISO_AT_FLIP.png');

fprintf('\n============================================================\n');
fprintf('STAGE 3 COMPLETE\n');
fprintf('Output folder:\n%s\n',outputFolder);
fprintf('Main questions answered by this run:\n');
fprintf('  1. Does SOC=0 also reorient at the same Ku?\n');
fprintf('  2. What is the Ku threshold in 0.25-0.50 meV?\n');
fprintf('  3. At what time does |Mz| cross 0.90?\n');
fprintf('  4. Which field dominates around that time?\n');
fprintf('Total wall time = %.2f min\n',toc(totalTimer)/60);
fprintf('============================================================\n');

end


%% =========================================================================
% NEW: PRECOMPUTE THE ANISOTROPY CASE (hoists the constant coefficient +
% unit vector out of the RHS -- this is the actual "efficiency" change)
% =========================================================================
function aniso = build_aniso_case_model(Ku_eV,easyAxis,surface_moment_eV_T)
aniso.Ku_eV = Ku_eV;
aniso.e_easy = easyAxis(:).'/norm(easyAxis);          % precomputed unit vector
aniso.coeff_T = 2*Ku_eV/surface_moment_eV_T;           % precomputed scalar coefficient
end


%% =========================================================================
% ONE CASE WITH ANISOTROPY: same two-stage integration as your control
% =========================================================================
function E = run_one_case_aniso(rho0,J0_meV,soc_meV,aniso,J_site_eV,exchange,model,bath,p,llgMode)
caseTimer=tic;
fprintf('\n------------------------------------------------------------\n');
fprintf('START J0=%g meV | SOC=%g meV | Ku=%g eV | %s\n',J0_meV,soc_meV,aniso.Ku_eV,llgMode);
fprintf('------------------------------------------------------------\n');

Nspin=size(model.H_static_spin,1);
Nrho=Nspin^2;
rhs=@(~,y) coupled_soc_exchange_aniso_rhs(y,J_site_eV,exchange,aniso,model,bath,p,llgMode);
y0=[rho0(:);p.M0(:)];

tObs=[]; MObs=[]; BmolObs=[]; BexObs=[]; BanisoObs=[]; BeffObs=[]; spinObs=[];
BmolPerpObs=[]; BexPerpObs=[]; BeffPerpObs=[]; spinPerpObs=[];
dMmagObs=[]; dMzObs=[]; tiltObs=[]; traceObs=[]; hermObs=[]; MnormObs=[];

% UNCHANGED: identical two-stage integration to your uploaded control.
optsEarly=odeset('InitialStep',p.InitialStep_fs,'MaxStep',p.MaxStepEarly_fs, ...
    'RelTol',p.RelTol,'AbsTol',p.AbsTol,'Refine',1,'Stats','off');
[tSeg,YSeg]=ode45(rhs,[p.tStart_fs p.earlyEnd_fs],y0,optsEarly);
[tObs,MObs,BmolObs,BexObs,BanisoObs,BeffObs,spinObs, ...
 BmolPerpObs,BexPerpObs,BeffPerpObs,spinPerpObs, ...
 dMmagObs,dMzObs,tiltObs,traceObs,hermObs,MnormObs,lastState] = ...
    append_full_segment_aniso(tSeg,YSeg,false,tObs,MObs,BmolObs,BexObs,BanisoObs,BeffObs,spinObs, ...
    BmolPerpObs,BexPerpObs,BeffPerpObs,spinPerpObs,dMmagObs,dMzObs,tiltObs, ...
    traceObs,hermObs,MnormObs,J_site_eV,aniso,model,p,llgMode);

optsLate=odeset('RelTol',p.RelTol,'AbsTol',p.AbsTol,'Refine',1,'Stats','off');
currentTime=p.earlyEnd_fs;
while currentTime<p.tEnd_fs-1e-12
    nextEnd=min(currentTime+p.chunk_fs,p.tEnd_fs);
    segTimer=tic;
    [tSeg,YSeg]=ode45(rhs,[currentTime nextEnd],lastState,optsLate);
    [tObs,MObs,BmolObs,BexObs,BanisoObs,BeffObs,spinObs, ...
     BmolPerpObs,BexPerpObs,BeffPerpObs,spinPerpObs, ...
     dMmagObs,dMzObs,tiltObs,traceObs,hermObs,MnormObs,lastState] = ...
        append_full_segment_aniso(tSeg,YSeg,true,tObs,MObs,BmolObs,BexObs,BanisoObs,BeffObs,spinObs, ...
        BmolPerpObs,BexPerpObs,BeffPerpObs,spinPerpObs,dMmagObs,dMzObs,tiltObs, ...
        traceObs,hermObs,MnormObs,J_site_eV,aniso,model,p,llgMode);
    currentTime=nextEnd;
    dt=diff(tSeg); if isempty(dt), dtMed=NaN; else, dtMed=median(dt); end
    fprintf('  Ku=%g: reached %5.2f ps | points=%7d | median dt=%g fs | %.2f min\n', ...
        aniso.Ku_eV,currentTime/1000,numel(tSeg),dtMed,toc(segTimer)/60);
end

rhoFinal=reshape(lastState(1:Nrho),Nspin,Nspin);
rhoFinal=0.5*(rhoFinal+rhoFinal');

E.J0_meV=J0_meV;
E.lambdaSOC_meV=soc_meV;
E.Ku_eV=aniso.Ku_eV;
E.t_fs=tObs;
E.M=MObs;
E.B_mol_mT=BmolObs;
E.B_ex_T=BexObs;
E.B_aniso_T=BanisoObs;
E.B_eff_T=BeffObs;
E.B_LLG_T=BeffObs;   % FULLXYZ: the full vector is what actually drives LLG
E.llgMode=llgMode;
E.totalElectronicSpin=spinObs;
E.Bmol_perp_T=BmolPerpObs;
E.Bex_perp_T=BexPerpObs;
E.Beff_perp_T=BeffPerpObs;
E.spin_perp=spinPerpObs;
E.dM_mag_per_fs=dMmagObs;
E.dMz_per_fs=dMzObs;
E.tilt_deg=tiltObs;
E.traceError=traceObs;
E.hermiticityError=hermObs;
E.MnormError=MnormObs;
E.finalRho=rhoFinal;
E.finalObs=density_observables_from_rho(rhoFinal,p);
E.runtime_min=toc(caseTimer)/60;

fprintf('COMPLETE Ku=%g | %.2f min | final M=[%.6g %.6g %.6g] | max|B_aniso|=%.6e T\n', ...
    aniso.Ku_eV,E.runtime_min,E.M(end,1),E.M(end,2),E.M(end,3), ...
    max(vecnorm(E.B_aniso_T,2,2)));
end


%% =========================================================================
% RHS WITH ANISOTROPY: only the B_eff sum changed vs. your control's RHS
% =========================================================================
function dy = coupled_soc_exchange_aniso_rhs(y,J_site_eV,exchange,aniso,model,bath,p,llgMode)
Nspin=size(model.H_static_spin,1);
Nrho=Nspin^2;
rho=reshape(y(1:Nrho),Nspin,Nspin);
rho=0.5*(rho+rho');
M=real(y(Nrho+1:Nrho+3)); M=M(:).'; M=M/max(norm(M),1e-15);

J_eV=compute_bond_current_general(rho,model);
B_mol_mT=apply_BS_kernel(model.BS_kernel_Rstar,J_eV);
B_mol_T=1e-3*B_mol_mT;

H_exchange=M(1)*exchange.Hx+M(2)*exchange.Hy+M(3)*exchange.Hz;
H_total=model.H_static_spin+H_exchange;

drho_H=-1i/p.hbar_eV_fs*(H_total*rho-rho*H_total);
drho_D=apply_fast_lindblad(rho,bath);
drho=drho_H+drho_D;

B_ex_T=compute_Bex_fast(rho,J_site_eV,model,p);

% Only new physics: anisotropy field, using the PRECOMPUTED coeff/axis
% (no norm() or division here -- both were done once in
% build_aniso_case_model, not on every RHS call).
B_aniso_T = aniso.coeff_T * dot(M,aniso.e_easy) * aniso.e_easy;

B_eff_T = B_mol_T + B_ex_T + B_aniso_T;

if strcmpi(llgMode,'ZONLY')
    B_LLG_T=[0 0 B_eff_T(3)];
elseif strcmpi(llgMode,'FULLXYZ')
    B_LLG_T=B_eff_T;
else
    error('Unknown llgMode: %s',llgMode);
end
precession=cross(M,B_LLG_T);
damping=cross(M,cross(M,B_LLG_T));
dM=p.gamma0_fs_T*(precession-p.lambda*damping);

dy=[drho(:);dM(:)];
end


%% =========================================================================
% APPEND EVERY ACCEPTED SOLVER POINT -- WITH B_aniso ALSO SAVED
% (this is exactly what your original append_full_segment did, extended
% with one extra field so the saved diagnostics actually reflect Ku)
% =========================================================================
function [tObs,MObs,BmolObs,BexObs,BanisoObs,BeffObs,spinObs, ...
          BmolPerpObs,BexPerpObs,BeffPerpObs,spinPerpObs, ...
          dMmagObs,dMzObs,tiltObs,traceObs,hermObs,MnormObs,lastState] = ...
    append_full_segment_aniso(tSeg,YSeg,dropFirst,tObs,MObs,BmolObs,BexObs,BanisoObs,BeffObs,spinObs, ...
    BmolPerpObs,BexPerpObs,BeffPerpObs,spinPerpObs,dMmagObs,dMzObs,tiltObs, ...
    traceObs,hermObs,MnormObs,J_site_eV,aniso,model,p,llgMode)

if dropFirst
    tSeg=tSeg(2:end); YSeg=YSeg(2:end,:);
end
n=numel(tSeg);
Nspin=size(model.H_static_spin,1); Nrho=Nspin^2;

Mseg=zeros(n,3); BmolSeg=zeros(n,3); BexSeg=zeros(n,3); BanisoSeg=zeros(n,3);
BeffSeg=zeros(n,3); spinSeg=zeros(n,3); BmolPerpSeg=zeros(n,1); BexPerpSeg=zeros(n,1);
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
    Baniso_T = aniso.coeff_T * dot(M,aniso.e_easy) * aniso.e_easy;
    Beff_T=Bmol_T+Bex_T+Baniso_T;

    obs=density_observables_from_rho(rho,p);
    sTot=[sum(obs.m_x),sum(obs.m_y),sum(obs.m_z)];

    BmolPerpVec=Bmol_T-dot(Bmol_T,M)*M;
    BexPerpVec=Bex_T-dot(Bex_T,M)*M;
    BeffPerpVec=Beff_T-dot(Beff_T,M)*M;
    spinPerpVec=sTot-dot(sTot,M)*M;

    if strcmpi(llgMode,'ZONLY')
        B_LLG_T=[0 0 Beff_T(3)];
    elseif strcmpi(llgMode,'FULLXYZ')
        B_LLG_T=Beff_T;
    else
        error('Unknown llgMode: %s',llgMode);
    end
    precession=cross(M,B_LLG_T);
    damping=cross(M,cross(M,B_LLG_T));
    dM=p.gamma0_fs_T*(precession-p.lambda*damping);

    Mseg(it,:)=M; BmolSeg(it,:)=Bmol_mT; BexSeg(it,:)=Bex_T; BanisoSeg(it,:)=Baniso_T;
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
BexObs=[BexObs;BexSeg]; BanisoObs=[BanisoObs;BanisoSeg]; BeffObs=[BeffObs;BeffSeg];
spinObs=[spinObs;spinSeg];
BmolPerpObs=[BmolPerpObs;BmolPerpSeg]; BexPerpObs=[BexPerpObs;BexPerpSeg];
BeffPerpObs=[BeffPerpObs;BeffPerpSeg]; spinPerpObs=[spinPerpObs;spinPerpSeg];
dMmagObs=[dMmagObs;dMmagSeg]; dMzObs=[dMzObs;dMzSeg]; tiltObs=[tiltObs;tiltSeg];
traceObs=[traceObs;traceSeg]; hermObs=[hermObs;hermSeg]; MnormObs=[MnormObs;MnormSeg];
lastState=YSeg(end,:).';
end


%% =========================================================================
% UNCHANGED HELPER FUNCTIONS -- copied verbatim from your uploaded script
% =========================================================================

function [Hunit,RhatVec,nnPairs] = build_board_nn_soc_unit(model)
L = size(model.Ri,1);
N = 2*L;
Nb = size(model.bonds,1);
sx=[0 1;1 0]; sy=[0 -1i;1i 0]; sz=[1 0;0 -1];
Sx=0.5*sx; Sy=0.5*sy; Sz=0.5*sz;
Hunit = complex(zeros(N,N));
RhatVec = zeros(Nb,3);
nnPairs = model.bonds;
for b = 1:Nb
    i = model.bonds(b,1); j = model.bonds(b,2);
    dR = model.Ri(j,:) - model.Ri(i,:);
    Rhat = dR/norm(dR);
    RhatVec(b,:) = Rhat;
    SdotR = Rhat(1)*Sx + Rhat(2)*Sy + Rhat(3)*Sz;
    block = 1i*SdotR;
    idxi = (2*i-1):(2*i); idxj = (2*j-1):(2*j);
    Hunit(idxi,idxj) = Hunit(idxi,idxj) + block;
    Hunit(idxj,idxi) = Hunit(idxj,idxi) + block';
end
end

function caseModel = build_soc_case_model(model,soc_meV)
caseModel = model;
lambdaSOC_eV = 1e-3*soc_meV;
caseModel.lambdaSOC_meV = soc_meV;
caseModel.H_static_spin = model.H_surface_spin_noSOC + lambdaSOC_eV*model.Hsoc_unit;
Nb = size(model.bonds,1);
caseModel.bondHij = complex(zeros(2,2,Nb));
for b = 1:Nb
    i = model.bonds(b,1); j = model.bonds(b,2);
    idxi = (2*i-1):(2*i); idxj = (2*j-1):(2*j);
    caseModel.bondHij(:,:,b) = caseModel.H_static_spin(idxi,idxj);
end
end

function J_eV = compute_bond_current_general(rho,model)
Nb = size(model.bonds,1);
J_eV = zeros(1,Nb);
for b = 1:Nb
    i = model.bonds(b,1); j = model.bonds(b,2);
    idxi = (2*i-1):(2*i); idxj = (2*j-1):(2*j);
    Hij = model.bondHij(:,:,b);
    rhoji = rho(idxj,idxi);
    J_eV(b) = -2*imag(trace(Hij*rhoji));
end
end

function J_eV = compute_bond_current_scalar_reference(rho,model)
Nb = size(model.bonds,1);
J_eV = zeros(1,Nb);
for b = 1:Nb
    i = model.bonds(b,1); j = model.bonds(b,2);
    tij = model.t_bonds_eV(b);
    up_i=2*i-1; dn_i=2*i; up_j=2*j-1; dn_j=2*j;
    J_eV(b) = 2*tij*imag(rho(up_j,up_i)+rho(dn_j,dn_i));
end
end

function J_site_eV = exchange_profile(J0_meV,model,p)
distance_A = vecnorm(model.Ri-p.Rstar_A,2,2);
J_site_eV = 1e-3*J0_meV*exp(-distance_A/p.r0_A);
end

function ex = precompute_exchange_matrices(J_site_eV)
L = numel(J_site_eV);
N = 2*L;
sx=[0 1;1 0]; sy=[0 -1i;1i 0]; sz=[1 0;0 -1];
Hx=complex(zeros(N,N)); Hy=complex(zeros(N,N)); Hz=complex(zeros(N,N));
for i = 1:L
    idx = (2*i-1):(2*i);
    Hx(idx,idx)=0.5*J_site_eV(i)*sx;
    Hy(idx,idx)=0.5*J_site_eV(i)*sy;
    Hz(idx,idx)=0.5*J_site_eV(i)*sz;
end
ex.Hx=Hx; ex.Hy=Hy; ex.Hz=Hz;
end

function B_ex_T = compute_Bex_fast(rho,J_site_eV,model,p)
ud = rho(model.idxUD);
mx = 2*real(ud); my = -2*imag(ud);
mz = real(rho(model.idxDiagUp))-real(rho(model.idxDiagDn));
J = J_site_eV(:);
weightedSpin_eV = [sum(J.*mx(:)),sum(J.*my(:)),sum(J.*mz(:))];
B_ex_T = -weightedSpin_eV/(2*p.surface_moment_eV_T);
end

function bath = build_fast_lindblad_bath(H0,p)
L = size(H0,1); N = 2*L;
[U,D] = eig(H0);
[E,idx] = sort(real(diag(D)));
U = U(:,idx);
T = kron(U,eye(2));
W = zeros(N,N);
for n = 1:L
    for m = n+1:L
        gap_eV = E(m)-E(n);
        if gap_eV <= 1e-12, continue; end
        gamma_down = p.Gamma0_per_fs;
        gamma_up = p.Gamma0_per_fs*exp(-p.beta_eV_inv*gap_eV);
        for spin = 1:2
            a_n = 2*n-2+spin; a_m = 2*m-2+spin;
            W(a_m,a_n) = gamma_up; W(a_n,a_m) = gamma_down;
        end
    end
end
outRate = sum(W,1).';
bath.T=T; bath.Tdag=T'; bath.W=W; bath.outRate=outRate;
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

function model = precompute_fast_indices(model)
L = size(model.Ri,1); N = 2*L;
up = (1:2:N).'; dn = (2:2:N).';
model.up=up; model.dn=dn;
model.idxDiagUp = sub2ind([N N],up,up);
model.idxDiagDn = sub2ind([N N],dn,dn);
model.idxUD = sub2ind([N N],up,dn);
end

function obs = density_observables_from_rho(rho,p)
N = size(rho,1); up = 1:2:N; down = 2:2:N;
obs.n_up = real(diag(rho(up,up))).';
obs.n_down = real(diag(rho(down,down))).';
obs.n_total = obs.n_up+obs.n_down;
rho_ud = diag(rho(up,down)).';
obs.m_x = 2*real(rho_ud); obs.m_y = -2*imag(rho_ud); obs.m_z = obs.n_up-obs.n_down;
obs.P_x = nan(size(obs.n_total)); obs.P_y = nan(size(obs.n_total)); obs.P_z = nan(size(obs.n_total));
valid = obs.n_total > p.density_threshold;
obs.P_x(valid) = obs.m_x(valid)./obs.n_total(valid);
obs.P_y(valid) = obs.m_y(valid)./obs.n_total(valid);
obs.P_z(valid) = obs.m_z(valid)./obs.n_total(valid);
end

function [rho0,diagnostic] = build_gas_ground_initial_state(H_hop)
[U,D] = eig(H_hop);
[E,idx] = sort(real(diag(D)));
U = U(:,idx); psi0 = U(:,1);
rhoOrbital = psi0*psi0';
rho0 = kron(rhoOrbital,0.5*eye(2));
diagnostic.energy_eV = E(1);
diagnostic.residual = norm(H_hop*psi0-E(1)*psi0);
diagnostic.traceError = abs(real(trace(rho0))-1);
diagnostic.hermiticityError = norm(rho0-rho0','fro');
up = 1:2:size(rho0,1); down = 2:2:size(rho0,1);
diagnostic.initialTotalMz = real(sum(diag(rho0(up,up)))-sum(diag(rho0(down,down))));
assert(diagnostic.residual < 1e-10,'Initial GS residual too large.');
assert(diagnostic.traceError < 1e-12,'Initial trace error too large.');
assert(diagnostic.hermiticityError < 1e-12,'Initial rho is not Hermitian.');
assert(abs(diagnostic.initialTotalMz) < 1e-12,'Initial molecular spin is not unpolarized.');
end

function model = build_upright_helix_model()
Nsites = 13; sitesPerTurn = 8; bondLength_A = 1.54;
dphi = 2*pi/sitesPerTurn;
axialRisePerSite_A = bondLength_A/sqrt(3);
transverseChord_A = bondLength_A*sqrt(2/3);
helixRadius_A = transverseChord_A/(2*sin(dphi/2));
siteIndex = (0:Nsites-1).'; phi = siteIndex*dphi;
Ri = [helixRadius_A*(cos(phi)-1), helixRadius_A*sin(phi), 1.0+axialRisePerSite_A*siteIndex];
bonds = [(1:Nsites-1).' (2:Nsites).'];
bondLengths_A = vecnorm(diff(Ri,1,1),2,2);
assert(max(abs(bondLengths_A-bondLength_A)) < 1e-12,'Helix bond lengths are not exactly 1.54 A.');
b1=Ri(2,:)-Ri(1,:); b2=Ri(3,:)-Ri(2,:); b3=Ri(4,:)-Ri(3,:);
chiralityTripleProduct = dot(b1,cross(b2,b3));
assert(abs(chiralityTripleProduct) > 1e-8,'Helix lost 3D chirality.');
model.name = 'helical_upright_13site_8perturn_1p5turn';
model.Ri=Ri; model.bonds=bonds; model.t_bonds=ones(size(bonds,1),1);
model.sitesPerTurn=sitesPerTurn; model.nTurns=(Nsites-1)/sitesPerTurn;
model.helixRadius_A=helixRadius_A; model.pitch_A=sitesPerTurn*axialRisePerSite_A;
model.chiralityTripleProduct=chiralityTripleProduct;
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
    i=bonds(b,1); j=bonds(b,2); tij=t_bonds_eV(b);
    H_hop(i,j)=-tij; H_hop(j,i)=-tij;
end
end

function kernel = precompute_BS_kernel_at_points(Ri,bonds,Rpoints_A,Nseg,min_distance_A)
nPoints = size(Rpoints_A,1); Nb = size(bonds,1);
kernel = zeros(3*nPoints,Nb);
for activeBond = 1:Nb
    unit_J_eV = zeros(1,Nb); unit_J_eV(activeBond) = 1;
    Bunit_mT = zeros(nPoints,3);
    for ip = 1:nPoints
        Bunit_mT(ip,:) = biot_savart_at_point(Ri,bonds,unit_J_eV,Rpoints_A(ip,:),0,Nseg,min_distance_A);
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
e_C = 1.602176634e-19; hbar_eV_s = 6.582119569e-16;
I_A = (e_C/hbar_eV_s)*J_eV;
mu0_over_4pi = 1e-7;
B_T = [0 0 0];
for b = 1:size(bonds,1)
    if abs(I_A(b)) < 1e-30, continue; end
    i=bonds(b,1); j=bonds(b,2);
    if excludedSite > 0 && (i==excludedSite || j==excludedSite), continue; end
    r1_A=Ri(i,:); r2_A=Ri(j,:); dl_A=(r2_A-r1_A)/Nseg;
    for s = 1:Nseg
        midpoint_A = r1_A+(s-0.5)*dl_A;
        rvec_A = Robs_A-midpoint_A; r_A = norm(rvec_A);
        if r_A < min_distance_A, continue; end
        dl_m = 1e-10*dl_A; rvec_m = 1e-10*rvec_A; r_m = norm(rvec_m);
        dB_T = mu0_over_4pi*I_A(b)*cross(dl_m,rvec_m)/(r_m^3);
        B_T = B_T+dB_T;
    end
end
B_mT = 1e3*B_T;
end

function save_figure(fig,outputFolder,fileName)
savefig(fig,fullfile(outputFolder,strrep(fileName,'.png','.fig')));
exportgraphics(fig,fullfile(outputFolder,fileName),'Resolution',300);
end

function downloadsDir = get_downloads_directory()
if ispc
    downloadsDir = fullfile(getenv('USERPROFILE'),'Downloads');
else
    homeDir = getenv('HOME');
    downloadsDir = fullfile(homeDir,'Downloads');
end
if ~exist(downloadsDir,'dir'), downloadsDir = pwd; end
end
