%% plot_hf14_zoom.m
% Loads the latest half-filling run and creates focused plots.
% Same data/units throughout; the By/Bz figure only changes axis limits.

downloadsDir = fullfile(getenv('USERPROFILE'),'Downloads');

files = dir(fullfile(downloadsDir,'**', ...
    'HALF_FILLING_N14_NO_SOC_J0_COMPARISON_COMPLETE.mat'));

assert(~isempty(files),'No HALF_FILLING_N14 result MAT file was found.');

[~,idx] = max([files.datenum]);
S = load(fullfile(files(idx).folder,files(idx).name),'runs');
runs = S.runs;

E0   = runs{find(cellfun(@(x) x.J0_meV==0,   runs),1)};
E150 = runs{find(cellfun(@(x) x.J0_meV==150, runs),1)};

B0   = 1e3*E0.B_eff_T;       % T -> mT
B150 = 1e3*E150.B_eff_T;     % T -> mT

% Focus windows (ps)
tZoomJ0   = 0.35;
tZoomJ150 = 1.50;
tZoomM    = 1.50;

%% FIGURE 1: J0=0 all effective-field components, tighter early-time zoom
figure('Color','w','Position',[180 100 1000 650]); hold on;
plot(E0.t_fs/1000,B0(:,1),'LineWidth',1.5,'DisplayName','B_x');
plot(E0.t_fs/1000,B0(:,2),'LineWidth',1.5,'DisplayName','B_y');
plot(E0.t_fs/1000,B0(:,3),'LineWidth',1.5,'DisplayName','B_z');
yline(0,'k--');
grid on;
xlabel('time [ps]');
ylabel('B_{eff} [mT]');
title('Half filling N=14 | J_0 = 0 meV | SOC=0');
legend('Location','best');
xlim([0 tZoomJ0]);

%% FIGURE 2: J0=150 all field components, zoom through settling
figure('Color','w','Position',[180 100 1000 650]); hold on;
plot(E150.t_fs/1000,B150(:,1),'LineWidth',1.5,'DisplayName','B_x');
plot(E150.t_fs/1000,B150(:,2),'LineWidth',1.5,'DisplayName','B_y');
plot(E150.t_fs/1000,B150(:,3),'LineWidth',1.5,'DisplayName','B_z');
yline(0,'k--');
grid on;
xlabel('time [ps]');
ylabel('B_{eff} [mT]');
title('Half filling N=14 | J_0 = 150 meV | SOC=0');
legend('Location','best');
xlim([0 tZoomJ150]);

%% FIGURE 3: J0=150 split scale (NO rescaling of data)
% Top: Bx in Tesla. Bottom: By and Bz in mT.
% This is only to make the much smaller Y/Z components visible.
figure('Color','w','Position',[180 70 1000 820]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile; hold on;
plot(E150.t_fs/1000,E150.B_eff_T(:,1),'LineWidth',1.6,'DisplayName','B_x');
yline(0,'k--');
grid on;
xlabel('time [ps]');
ylabel('B_x [T]');
title('J_0 = 150 meV | dominant longitudinal component');
legend('Location','best');
xlim([0 tZoomJ150]);

nexttile; hold on;
plot(E150.t_fs/1000,B150(:,2),'LineWidth',1.7,'DisplayName','B_y');
plot(E150.t_fs/1000,B150(:,3),'LineWidth',1.7,'DisplayName','B_z');
yline(0,'k--');
grid on;
xlabel('time [ps]');
ylabel('B_{y,z} [mT]');
title('Same J_0 = 150 run | zoom on smaller Y/Z components');
legend('Location','best');
xlim([0 tZoomJ150]);

%% FIGURE 4: Mx / My / Mz separately, both J0 models, early-time zoom
figure('Color','w','Position',[130 60 1400 900]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

labels = {'M_x','M_y','M_z'};

for comp = 1:3
    nexttile; hold on;

    plot(E0.t_fs/1000,E0.M(:,comp),'LineWidth',1.55, ...
        'DisplayName','J_0 = 0');
    plot(E150.t_fs/1000,E150.M(:,comp),'LineWidth',1.55, ...
        'DisplayName','J_0 = 150 meV');

    yline(0,'k--');
    grid on;
    xlim([0 tZoomM]);
    ylabel(labels{comp});
    title(labels{comp});
    legend('Location','best');

    % Each component gets its own y-axis autoscale.
    idx0   = E0.t_fs/1000   <= tZoomM;
    idx150 = E150.t_fs/1000 <= tZoomM;
    vals = [E0.M(idx0,comp); E150.M(idx150,comp)];

    ymin = min(vals);
    ymax = max(vals);
    span = ymax-ymin;

    if span > 1e-14
        pad = 0.08*span;
        ylim([ymin-pad ymax+pad]);
    end
end

xlabel('time [ps]');
sgtitle('Half filling N=14 | SOC=0 | early-time moment response');
