function no_soc_xyz_control()
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

% CLEAN CONTROL: exchange stays ON exactly as in the current J0=150 model.
p.J0List_meV = [150];
p.r0_A = 3.0;

% CLEAN CONTROL: SOC is the ONLY Hamiltonian ingredient switched off.
p.lambdaSOCList_meV = [0];

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
p.requestedWorkers = 1;

assert(isequal(p.M0,[1 0 0]),'M(0) must be +X.');
assert(abs(p.lambda-0.6)<1e-14,'lambda_LLG must remain 0.6.');
assert(abs(p.tau_rel_fs-500)<1e-12,'tau_rel must remain 500 fs.');
assert(any(p.lambdaSOCList_meV==0),'SOC=0 control is required.');

fprintf('\n============================================================\n');
fprintf('NO-SOC CONTROL | J0=150 meV | Z-ONLY vs FULL XYZ LLG\n');
fprintf('============================================================\n');
fprintf('J0 scan           = '); fprintf('%g ',p.J0List_meV); fprintf('meV\n');
fprintf('SOC scan          = '); fprintf('%g ',p.lambdaSOCList_meV); fprintf('meV\n');
fprintf('M(0)              = [%g %g %g]\n',p.M0);
fprintf('lambda_LLG        = %.3f\n',p.lambda);
fprintf('tau_rel           = %.1f fs\n',p.tau_rel_fs);
fprintf('physical time     = %.2f ps\n',p.tEnd_fs/1000);
fprintf('SOC               = 0 meV exactly (H_SOC term = 0)\n');
fprintf('LLG fields tested = Z-only and FULL XYZ\n');
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
    ['NO_SOC_EXCHANGE_J150_ZONLY_VS_FULLXYZ_10PS_' runTag]);
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


%% 7. TWO CLEAN NO-SOC RUNS =============================================
% Everything below is identical between the two runs except the vector
% supplied to LLG:
%
%   ZONLY   : B_LLG = [0, 0, B_eff,z]
%   FULLXYZ : B_LLG = [B_eff,x, B_eff,y, B_eff,z]
%
% In BOTH:
%   J0 = 150 meV
%   SOC = 0 meV
%   M(0) = +X
%   lambda_LLG = 0.6
%   tau_rel = 500 fs
%   runtime = 10 ps
%   same geometry, rho0, bath, R*, exchange profile, solver tolerances
%   electronic Zeeman = OFF

J0_meV = 150;
soc_meV = 0;

assert(numel(p.J0List_meV)==1 && p.J0List_meV==150, ...
    'This control must keep J0=150 meV.');
assert(numel(p.lambdaSOCList_meV)==1 && p.lambdaSOCList_meV==0, ...
    'This control must have SOC exactly zero.');

J_site_eV = JsiteCell{1};
exchange = exchangeCell{1};

E_Z = run_one_case(rho0,J0_meV,soc_meV,J_site_eV,exchange, ...
    model,bath,p,'ZONLY');

E_XYZ = run_one_case(rho0,J0_meV,soc_meV,J_site_eV,exchange, ...
    model,bath,p,'FULLXYZ');

%% 8. AUDIT: CONFIRM ONLY LLG FIELD SELECTION CHANGED ====================
fprintf('\n============================================================\n');
fprintf('CONTROL AUDIT\n');
fprintf('============================================================\n');
fprintf('J0                 : %.1f meV in BOTH\n',J0_meV);
fprintf('SOC                : %.1f meV in BOTH\n',soc_meV);
fprintf('M0                 : [%g %g %g] in BOTH\n',p.M0);
fprintf('lambda_LLG         : %.3f in BOTH\n',p.lambda);
fprintf('tau_rel            : %.1f fs in BOTH\n',p.tau_rel_fs);
fprintf('runtime            : %.2f ps in BOTH\n',p.tEnd_fs/1000);
fprintf('R*                 : [%.6f %.6f %.6f] A in BOTH\n',p.Rstar_A);
fprintf('electronic Zeeman  : OFF in BOTH\n');
fprintf('ONLY difference    : B_LLG selection\n');
fprintf('  ZONLY             B_LLG=(0,0,B_eff,z)\n');
fprintf('  FULLXYZ           B_LLG=(B_eff,x,B_eff,y,B_eff,z)\n');
fprintf('============================================================\n');

%% 9. FIGURE A: ACTUAL LLG FIELD =========================================
figA = figure('Color','w','Position',[220 70 1550 820]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile; hold on;
plot(E_Z.t_fs/1000,1e3*E_Z.B_LLG_T(:,1),'LineWidth',1.2,'DisplayName','B_{LLG,x}');
plot(E_Z.t_fs/1000,1e3*E_Z.B_LLG_T(:,2),'LineWidth',1.2,'DisplayName','B_{LLG,y}');
plot(E_Z.t_fs/1000,1e3*E_Z.B_LLG_T(:,3),'LineWidth',1.2,'DisplayName','B_{LLG,z}');
grid on;
xlabel('time [ps]'); ylabel('B_{LLG} [mT]');
title('NO SOC | J_0=150 meV | Z-only LLG');
legend('Location','best');

nexttile; hold on;
plot(E_XYZ.t_fs/1000,1e3*E_XYZ.B_LLG_T(:,1),'LineWidth',1.2,'DisplayName','B_{LLG,x}');
plot(E_XYZ.t_fs/1000,1e3*E_XYZ.B_LLG_T(:,2),'LineWidth',1.2,'DisplayName','B_{LLG,y}');
plot(E_XYZ.t_fs/1000,1e3*E_XYZ.B_LLG_T(:,3),'LineWidth',1.2,'DisplayName','B_{LLG,z}');
grid on;
xlabel('time [ps]'); ylabel('B_{LLG} [mT]');
title('NO SOC | J_0=150 meV | Full XYZ LLG');
legend('Location','best');

sgtitle('Actual field entering LLG — identical model, SOC = 0');
save_figure(figA,outputFolder,'FIG_NO_SOC_J150_FIELD_ZONLY_VS_FULLXYZ.png');

%% 10. FIGURE B: Mz DIRECT COMPARISON ====================================
figB = figure('Color','w','Position',[260 110 1250 720]); hold on;
plot(E_Z.t_fs/1000,E_Z.M(:,3),'LineWidth',1.55, ...
    'DisplayName','SOC=0 | Z-only LLG');
plot(E_XYZ.t_fs/1000,E_XYZ.M(:,3),'LineWidth',1.55, ...
    'DisplayName','SOC=0 | Full XYZ LLG');
yline(0,'--');
grid on;
xlabel('time [ps]');
ylabel('M_z');
title('NO SOC | J_0=150 meV | M_z: Z-only vs Full XYZ');
legend('Location','best');
xlim([0 p.tEnd_fs/1000]);
save_figure(figB,outputFolder,'FIG_NO_SOC_J150_MZ_ZONLY_VS_FULLXYZ.png');

%% 11. FIGURE C: PHYSICAL Beff COMPONENTS ================================
% This is useful because Z-only still has physical B_eff,x/y; they are
% simply excluded from the LLG torque.  Plotting B_eff separately prevents
% confusion between "field exists" and "field enters LLG".
figC = figure('Color','w','Position',[220 70 1550 820]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile; hold on;
plot(E_Z.t_fs/1000,1e3*E_Z.B_eff_T(:,1),'LineWidth',1.2,'DisplayName','B_{eff,x}');
plot(E_Z.t_fs/1000,1e3*E_Z.B_eff_T(:,2),'LineWidth',1.2,'DisplayName','B_{eff,y}');
plot(E_Z.t_fs/1000,1e3*E_Z.B_eff_T(:,3),'LineWidth',1.2,'DisplayName','B_{eff,z}');
grid on;
xlabel('time [ps]'); ylabel('B_{eff} [mT]');
title('SOC=0 | trajectory generated with Z-only LLG');
legend('Location','best');

nexttile; hold on;
plot(E_XYZ.t_fs/1000,1e3*E_XYZ.B_eff_T(:,1),'LineWidth',1.2,'DisplayName','B_{eff,x}');
plot(E_XYZ.t_fs/1000,1e3*E_XYZ.B_eff_T(:,2),'LineWidth',1.2,'DisplayName','B_{eff,y}');
plot(E_XYZ.t_fs/1000,1e3*E_XYZ.B_eff_T(:,3),'LineWidth',1.2,'DisplayName','B_{eff,z}');
grid on;
xlabel('time [ps]'); ylabel('B_{eff} [mT]');
title('SOC=0 | trajectory generated with Full XYZ LLG');
legend('Location','best');

sgtitle('Physical effective field B_{eff}=B_{mol}+B_{ex} — SOC = 0');
save_figure(figC,outputFolder,'FIG_NO_SOC_J150_BEFF_COMPONENTS_ZONLY_VS_FULLXYZ.png');

%% 12. NUMERIC SUMMARY ====================================================
mode = {'Z-only';'Full XYZ'};
Mz_end = [E_Z.M(end,3); E_XYZ.M(end,3)];
max_abs_Mz = [max(abs(E_Z.M(:,3))); max(abs(E_XYZ.M(:,3)))];
max_Beff_perp_mT = 1e3*[max(E_Z.Beff_perp_T); max(E_XYZ.Beff_perp_T)];
max_abs_Beff_z_mT = 1e3*[max(abs(E_Z.B_eff_T(:,3))); ...
                         max(abs(E_XYZ.B_eff_T(:,3)))];
max_abs_BLLG_mT = 1e3*[max(vecnorm(E_Z.B_LLG_T,2,2)); ...
                       max(vecnorm(E_XYZ.B_LLG_T,2,2))];

summaryTable = table(mode, ...
    repmat(J0_meV,2,1),repmat(soc_meV,2,1), ...
    Mz_end,max_abs_Mz,max_Beff_perp_mT,max_abs_Beff_z_mT,max_abs_BLLG_mT, ...
    'VariableNames',{'LLG_mode','J0_meV','SOC_meV','Mz_end','max_abs_Mz', ...
    'max_Beff_perp_mT','max_abs_Beff_z_mT','max_abs_BLLG_mT'});

disp(summaryTable);
writetable(summaryTable,fullfile(outputFolder, ...
    'NO_SOC_J150_ZONLY_VS_FULLXYZ_SUMMARY.csv'));

save(fullfile(outputFolder,'NO_SOC_J150_ZONLY_VS_FULLXYZ_COMPLETE.mat'), ...
    'E_Z','E_XYZ','summaryTable','p','model','initialDiagnostic','-v7.3');

fprintf('\nSaved clean no-SOC comparison to:\n%s\n',outputFolder);
fprintf('Total runtime: %.2f min\n',toc(totalTimer)/60);

end

%% =========================================================================
% ONE CASE: Z-ONLY OR FULL XYZ
% =========================================================================
function E = run_one_case(rho0,J0_meV,soc_meV,J_site_eV,exchange,model,bath,p,llgMode)
caseTimer=tic;
fprintf('\n------------------------------------------------------------\n');
fprintf('START J0=%g meV | SOC=%g meV | %s\n',J0_meV,soc_meV,llgMode);
fprintf('------------------------------------------------------------\n');

caseModel=build_soc_case_model(model,soc_meV);
Nspin=size(caseModel.H_static_spin,1);
Nrho=Nspin^2;
rhs=@(~,y) coupled_soc_exchange_rhs_mode(y,J_site_eV,exchange,caseModel,bath,p,llgMode);
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
    traceObs,hermObs,MnormObs,J_site_eV,caseModel,p,llgMode);

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
        traceObs,hermObs,MnormObs,J_site_eV,caseModel,p,llgMode);
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
if strcmpi(llgMode,'ZONLY')
    E.B_LLG_T=[zeros(size(BeffObs,1),2),BeffObs(:,3)];
elseif strcmpi(llgMode,'FULLXYZ')
    E.B_LLG_T=BeffObs;
else
    error('Unknown llgMode: %s',llgMode);
end
E.llgMode=llgMode;
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

fprintf('COMPLETE %s | J0=%g SOC=%g | %.2f min | final M=[%.6g %.6g %.6g] | max Beff_perp=%.6e T\n', ...
    llgMode,J0_meV,soc_meV,E.runtime_min,E.M(end,1),E.M(end,2),E.M(end,3),max(E.Beff_perp_T));
end


%% =========================================================================
% COUPLED RHS: SELECTED LLG FIELD
% =========================================================================
function dy = coupled_soc_exchange_rhs_mode(y,J_site_eV,exchange,model,bath,p,llgMode)
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

% The ONLY controlled difference between the two no-SOC runs.
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
% APPEND EVERY ACCEPTED SOLVER POINT
% =========================================================================
function [tObs,MObs,BmolObs,BexObs,BeffObs,spinObs, ...
          BmolPerpObs,BexPerpObs,BeffPerpObs,spinPerpObs, ...
          dMmagObs,dMzObs,tiltObs,traceObs,hermObs,MnormObs,lastState] = ...
    append_full_segment(tSeg,YSeg,dropFirst,tObs,MObs,BmolObs,BexObs,BeffObs,spinObs, ...
    BmolPerpObs,BexPerpObs,BeffPerpObs,spinPerpObs,dMmagObs,dMzObs,tiltObs, ...
    traceObs,hermObs,MnormObs,J_site_eV,model,p,llgMode)

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

    % Use the exact same LLG field selection as in the ODE.
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
