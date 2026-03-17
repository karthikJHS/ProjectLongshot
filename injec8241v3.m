%% =========================================================
%  N2O COLD FLOW ANALYSIS TOOL  —  coldflow_analysis.m
%  =========================================================
%  Characterises nitrous oxide flow through your injector
%  system using load cell data from a cold flow test.
%  N2O only. No combustion. No fuel.
%
%  WHAT THIS SCRIPT DOES:
%    1. Loads your load cell CSV and plots the raw data
%    2. Computes N2O mass flow rate (dm/dt)
%    3. Branches into two analysis paths:
%
%       PATH A — Injector Characterisation
%         You have a specific drill bit size installed.
%         Evaluates how that orifice actually performed:
%         effective Cd, measured vs predicted flow, SPI vs Dyer.
%
%       PATH B — Injector Sizing
%         You want to find the best drill bit size.
%         Sweeps a range of standard drill bit diameters and
%         recommends the best match to your measured N2O flow.
%
%  CSV FORMAT (minimum):
%    Column 1: time          [s]
%    Column 2: mass          [kg or g]  — you will be asked
%
%  CSV FORMAT (with live sensor data, optional):
%    Column 1: time          [s]
%    Column 2: mass          [kg or g]
%    Column 3: P_inlet       [psi]  upstream of injector
%    Column 4: P_downstream  [psi]  downstream of injector
%    Column 5: D_inj         [in]   orifice diameter
%
%  SYSTEM DEFAULTS (from your pressure budget):
%    Tank:        627 - 653 psi  (15 deg C saturated N2O)
%    P_inlet:     577 psi        (after hose + needle valve losses)
%    P_chamber:   440 psi        (target hot fire chamber pressure)
%    Current inj: 18018NOS, D = 0.094 in, Cd = 0.7
% =========================================================

clear; clc; close all;

fprintf('╔══════════════════════════════════════════════╗\n');
fprintf('║   N2O COLD FLOW ANALYSIS TOOL                ║\n');
fprintf('║   Injector Characterisation & Sizing         ║\n');
fprintf('╚══════════════════════════════════════════════╝\n\n');

%% ── STEP 1: LOAD CSV ────────────────────────────────────
[t, mass_kg, csv_extras, filename] = load_csv_data();

%% ── STEP 2: PROCESS & PLOT RAW DATA ─────────────────────
[t_trim, ~, mdot, mdot_avg, mdot_peak] = process_and_plot(t, mass_kg, filename);

%% ── STEP 3: ANALYSIS MENU ───────────────────────────────
fprintf('\n╔══════════════════════════════════════════════╗\n');
fprintf('║   SELECT ANALYSIS PATH                       ║\n');
fprintf('╠══════════════════════════════════════════════╣\n');
fprintf('║  A  →  Injector Characterisation             ║\n');
fprintf('║         You have a drill bit size installed. ║\n');
fprintf('║         Evaluate how it actually performed.  ║\n');
fprintf('║                                              ║\n');
fprintf('║  B  →  Injector Sizing                       ║\n');
fprintf('║         Find the best drill bit size from    ║\n');
fprintf('║         your measured N2O flow.              ║\n');
fprintf('╚══════════════════════════════════════════════╝\n');

while true
    choice = upper(strtrim(input('\nEnter A or B: ', 's')));
    if ismember(choice, {'A','B'}), break; end
    fprintf('  Invalid input. Please enter A or B.\n');
end

%% ── STEP 4: BRANCH ───────────────────────────────────────
switch choice
    case 'A'
        fprintf('\n── PATH A: INJECTOR CHARACTERISATION ──────────\n');
        inj = get_injector_geometry(csv_extras);
        injector_characterisation(t_trim, mdot, mdot_avg, mdot_peak, inj);

    case 'B'
        fprintf('\n── PATH B: INJECTOR SIZING ─────────────────────\n');
        targets = get_target_conditions();
        injector_sizing(mdot_avg, mdot_peak, targets);
end

fprintf('\n✓ Cold flow analysis complete.\n\n');


%% =========================================================
%  LOCAL FUNCTIONS
%% =========================================================

% ---------------------------------------------------------
function [t, mass_kg, csv_extras, filename] = load_csv_data()
% Load CSV file, handle unit conversion, detect optional sensor columns.

fprintf('\n── STEP 1: LOAD CSV DATA ───────────────────────\n');

[file, path] = uigetfile('*.csv', 'Select Cold Flow CSV File');
if isequal(file, 0), error('No file selected. Exiting.'); end

filename = file;
raw      = readmatrix(fullfile(path, file));
fprintf('  Loaded: %s\n', filename);

if size(raw, 2) < 2
    error('CSV must have at least 2 columns: time and mass.');
end

t    = raw(:,1);
mass = raw(:,2);

fprintf('\n  What units is the mass column in?\n');
fprintf('    1) Kilograms [kg]\n');
fprintf('    2) Grams [g]\n');
fprintf('    3) Pounds [lb]\n');
while true
    u = strtrim(input('  Enter 1, 2, or 3: ', 's'));
    if ismember(u, {'1','2','3'}), break; end
end
switch u
    case '1', mass_kg = mass;
    case '2', mass_kg = mass / 1000;
    case '3', mass_kg = mass * 0.453592;
end

fprintf('  Time:  %.3f s  to  %.3f s\n', t(1), t(end));
fprintf('  Mass:  %.4f kg  to  %.4f kg\n', max(mass_kg), min(mass_kg));

% Detect optional extra sensor columns
csv_extras = struct();
if size(raw, 2) >= 3
    fprintf('\n  %d columns detected.\n', size(raw,2));
    fprintf('  Do columns 3-5 contain live sensor data?\n');
    fprintf('  [P_inlet (psi), P_downstream (psi), D_inj (in)]\n');
    if strcmp(upper(strtrim(input('  Y / N: ', 's'))), 'Y')
        if size(raw,2) >= 3, csv_extras.P_inlet      = raw(:,3); end
        if size(raw,2) >= 4, csv_extras.P_downstream = raw(:,4); end
        if size(raw,2) >= 5, csv_extras.D_inj        = raw(:,5); end
        fprintf('  Sensor columns loaded.\n');
    end
end
end

% ---------------------------------------------------------
function [t_trim, mass_trim, mdot, mdot_avg, mdot_peak] = process_and_plot(t, mass_kg, filename)
% Smooth mass signal, compute -dm/dt, auto-detect flow window, plot.

fprintf('\n── STEP 2: PROCESS & PLOT ──────────────────────\n');

% Smooth to reduce load cell noise before differentiating
smooth_window = max(5, round(length(t) * 0.02));
mass_smooth   = smoothdata(mass_kg, 'gaussian', smooth_window);

% Mass flow = -dm/dt (tank mass decreases as N2O flows out)
mdot_raw = -gradient(mass_smooth, t);

% Auto-detect flow window: where mdot exceeds 5% of its peak value
threshold    = 0.05 * max(mdot_raw);
flow_indices = find(mdot_raw > threshold);

if isempty(flow_indices)
    warning('Could not auto-detect flow window. Using full dataset.');
    flow_indices = 1:length(t);
end

i_start = flow_indices(1);
i_end   = flow_indices(end);

fprintf('  Auto-detected flow window: %.3f s  to  %.3f s\n', t(i_start), t(i_end));
fprintf('  Accept? (N to enter manually)\n');
if strcmp(upper(strtrim(input('  Y / N: ', 's'))), 'N')
    i_start = find(t >= input('  Start time [s]: '), 1, 'first');
    i_end   = find(t <= input('  End time   [s]: '), 1, 'last');
    fprintf('  Using: %.3f s  to  %.3f s\n', t(i_start), t(i_end));
end

t_trim    = t(i_start:i_end);
mass_trim = mass_smooth(i_start:i_end);
mdot      = mdot_raw(i_start:i_end);
mdot_avg  = mean(mdot);
mdot_peak = max(mdot);
total_n2o = mass_kg(i_start) - mass_kg(i_end);

fprintf('\n  ┌──────────────────────────────────────────────┐\n');
fprintf('  │  N2O MASS FLOW SUMMARY                       │\n');
fprintf('  ├──────────────────────────────────────────────┤\n');
fprintf('  │  Average mdot  :  %7.4f kg/s  (%6.2f g/s)  │\n', mdot_avg,  mdot_avg*1000);
fprintf('  │  Peak    mdot  :  %7.4f kg/s  (%6.2f g/s)  │\n', mdot_peak, mdot_peak*1000);
fprintf('  │  Total N2O used:  %7.4f kg    (%6.2f g)    │\n', total_n2o, total_n2o*1000);
fprintf('  │  Flow duration :  %7.3f s                   │\n', t_trim(end)-t_trim(1));
fprintf('  └──────────────────────────────────────────────┘\n');

figure('Name','Cold Flow  --  Raw Data','Position',[80 80 980 600]);

subplot(2,1,1);
plot(t, mass_kg*1000, 'Color',[0.75 0.75 0.75], 'LineWidth',1.2, 'DisplayName','Raw');
hold on;
plot(t, mass_smooth*1000, 'b-', 'LineWidth',2, 'DisplayName','Smoothed');
xline(t(i_start), '--k', 'LineWidth',1.3, 'DisplayName','Flow window');
xline(t(i_end),   '--k', 'LineWidth',1.3, 'HandleVisibility','off');
xlabel('Time [s]'); ylabel('N2O Mass [g]');
title(['Mass vs Time  --  ', strrep(filename,'_','\_')], 'FontSize',12);
legend('Location','northeast','FontSize',9); grid on;

subplot(2,1,2);
plot(t_trim, mdot*1000, 'r-', 'LineWidth',2, 'DisplayName','Measured mdot');
hold on;
yline(mdot_avg*1000,  '--b', 'LineWidth',1.8, 'DisplayName',sprintf('Avg:  %.2f g/s', mdot_avg*1000));
yline(mdot_peak*1000, ':m',  'LineWidth',1.8, 'DisplayName',sprintf('Peak: %.2f g/s', mdot_peak*1000));
xlabel('Time [s]'); ylabel('Mass Flow Rate [g/s]');
title('N2O Mass Flow Rate  (dm/dt)', 'FontSize',12);
legend('Location','northeast','FontSize',9); grid on;

sgtitle('Cold Flow Test  --  N2O Load Cell Data', 'FontSize',14, 'FontWeight','bold');
end

% ---------------------------------------------------------
function inj = get_injector_geometry(csv_extras)
% Collect injector orifice geometry: manual entry or from CSV.
% Accepts multiple drill bit sizes for side-by-side comparison.

fprintf('\n── INJECTOR GEOMETRY ───────────────────────────\n');

has_csv = isfield(csv_extras,'D_inj') && ~isempty(csv_extras.D_inj);

if has_csv
    fprintf('  Orifice diameter data found in CSV.\n');
    fprintf('    1) Use diameter from CSV\n');
    fprintf('    2) Enter manually\n');
    while true
        src = strtrim(input('  Enter 1 or 2: ','s'));
        if ismember(src,{'1','2'}), break; end
    end
else
    src = '2';
end

if strcmp(src,'1')
    inj.D_inj        = csv_extras.D_inj;
    inj.P_inlet      = csv_extras.P_inlet;
    inj.P_downstream = csv_extras.P_downstream;
    inj.source       = 'csv';
    fprintf('  D_inj mean: %.4f in\n', mean(inj.D_inj));
    fprintf('  P_inlet avg: %.1f psi\n', mean(inj.P_inlet));
    fprintf('  P_downstream avg: %.1f psi\n', mean(inj.P_downstream));
else
    fprintf('\n  Enter the drill bit size(s) you tested [inches].\n');
    fprintf('  Single:    0.094\n');
    fprintf('  Multiple:  [0.070, 0.094, 0.110]\n');
    fprintf('\n  Common drill bit reference:\n');
    fprintf('    #55=0.052  #52=0.063  #50=0.070  #48=0.076\n');
    fprintf('    #44=0.086  #42=0.094  #40=0.098  #38=0.104\n');
    fprintf('    #35=0.110  #31=0.120  #29=0.136  #26=0.147\n');
    inj.D_inj = input('\n  D_inj [in] = ');

    fprintf('\n  P_inlet [psi]  (pressure budget default: 577 psi)\n');
    inj.P_inlet = parse_default(strtrim(input('  P_inlet: ','s')), 577);

    fprintf('\n  P_downstream [psi]\n');
    fprintf('  Cold flow no chamber = 14.7 (atmospheric)\n');
    fprintf('  Chamber stand-in = enter that pressure\n');
    inj.P_downstream = parse_default(strtrim(input('  P_downstream: ','s')), 14.7);

    inj.source = 'manual';
end

inj.Cd         = parse_default(strtrim(input('\n  Discharge coeff Cd  (default 0.7): ','s')), 0.7);
inj.n_orifices = parse_default(strtrim(input('  Number of orifices  (default 1):   ','s')), 1);

D_m         = inj.D_inj * 0.0254;
inj.A_total = (pi/4 .* D_m.^2) * inj.n_orifices;

fprintf('\n  Injector summary:\n');
for k = 1:length(inj.D_inj)
    fprintf('    D=%.4f in  |  A_total=%.6f cm2  |  x%d orifices\n', ...
        inj.D_inj(k), inj.A_total(k)*1e4, inj.n_orifices);
end
fprintf('    Cd=%.2f  |  P_inlet=%.1f psi  |  P_downstream=%.1f psi\n', ...
    inj.Cd, mean(inj.P_inlet), mean(inj.P_downstream));
end

% ---------------------------------------------------------
function injector_characterisation(t_trim, mdot, mdot_avg, mdot_peak, inj)
% PATH A: Evaluate how a specific drill bit orifice performed.
%
% SPI  — Single-Phase Incompressible: assumes N2O stays liquid.
%         Simpler, slightly overestimates flow.
%
% Dyer — Two-phase blended model: accounts for N2O partially
%         flashing to vapour as it passes through the orifice.
%         More accurate for saturated N2O systems.
%
% Effective Cd: back-calculated from your actual measured flow.
%   - Cd_eff < design Cd  ==>  orifice is more restrictive than ideal
%   - Cd_eff > design Cd  ==>  orifice is less restrictive than ideal
%   - Cd_eff close to 0.7 ==>  your drilled hole is clean and sharp

fprintf('\n── PATH A: INJECTOR CHARACTERISATION ──────────\n');

% N2O fluid properties at 15 deg C saturated liquid
rho_l = 793;        % [kg/m3] liquid density
rho_g = 128;        % [kg/m3] vapour density
h_fg  = 323.6e3;    % [J/kg]  latent heat of vaporisation
P_sat = 5169e3;     % [Pa]    saturation pressure at 15 deg C

P_inlet_pa = mean(inj.P_inlet)      * 6894.76;
% [] Pressure from psi to pa
P_down_pa  = mean(inj.P_downstream) * 6894.76;
% [] Pressure from psi to pa
dP         = max(P_inlet_pa - P_down_pa, 0);
Cd         = inj.Cd;

% SPI prediction
%mdot_SPI = Cd .* inj.A_total .* sqrt(2 * rho_l * dP);
dh = h1 - h2;
V2 = sqrt(dh*2);
mdot = Cd * A * rho2 * V2;

% Dyer prediction
x_exit   = min(max(dP / (rho_l * h_fg), 0), 1);
rho_2ph  = 1 / (x_exit/rho_g + (1-x_exit)/rho_l);
v_hem    = sqrt(2 * dP / rho_l);
mdot_HEM = Cd .* inj.A_total .* rho_2ph .* v_hem;
if P_inlet_pa >= P_sat, kappa = 0;
else, kappa = sqrt(min(max(dP/max(P_inlet_pa-P_sat,1),0),1)); end
mdot_Dyer = (1-kappa).*mdot_SPI + kappa.*mdot_HEM;

% Effective Cd from measured flow
Cd_eff = mdot_avg ./ (inj.A_total .* sqrt(2 * rho_l * dP));

fprintf('\n  Pressure drop across orifice: %.2f psi\n\n', dP/6894.76);
fprintf('  ┌──────────┬──────────┬──────────┬──────────┬──────────┬───────────┐\n');
fprintf('  │ D [in]   │ Measured │ SPI pred │ Dyer pred│ Cd_eff   │ Dyer err  │\n');
fprintf('  │          │ avg [g/s]│ [g/s]    │ [g/s]    │          │ vs meas   │\n');
fprintf('  ├──────────┼──────────┼──────────┼──────────┼──────────┼───────────┤\n');
for k = 1:length(inj.D_inj)
    err = (mdot_Dyer(k) - mdot_avg) / mdot_avg * 100;
    fprintf('  │ D=%.4f │ %6.3f   │ %6.3f   │ %6.3f   │  %.4f  │  %+.1f%%   │\n', ...
        inj.D_inj(k), mdot_avg*1000, mdot_SPI(k)*1000, mdot_Dyer(k)*1000, Cd_eff(k), err);
end
fprintf('  └──────────┴──────────┴──────────┴──────────┴──────────┴───────────┘\n\n');

% Downstream pressure sweep: shows predicted N2O flow across
% full range from atmospheric (cold flow) to 440 psi (hot fire)
P_sweep_pa  = linspace(14.7*6894.76, 500*6894.76, 300);
P_sweep_psi = P_sweep_pa / 6894.76;
n_d = length(inj.D_inj);
mdot_SPI_sw  = zeros(n_d, length(P_sweep_pa));
mdot_Dyer_sw = zeros(n_d, length(P_sweep_pa));

for k = 1:n_d
    for i = 1:length(P_sweep_pa)
        dPi   = max(P_inlet_pa - P_sweep_pa(i), 0);
        spi_i = Cd * inj.A_total(k) * sqrt(2*rho_l*dPi);
        mdot_SPI_sw(k,i) = spi_i;
        x_i   = min(max(dPi/(rho_l*h_fg),0),1);
        rho_i = 1/(x_i/rho_g + (1-x_i)/rho_l);
        hem_i = Cd * inj.A_total(k) * rho_i * sqrt(2*dPi/rho_l);
        if P_inlet_pa >= P_sat, ki=0;
        else, ki=sqrt(min(max(dPi/max(P_inlet_pa-P_sat,1),0),1)); end
        mdot_Dyer_sw(k,i) = (1-ki)*spi_i + ki*hem_i;
    end
end

colors = lines(n_d);
figure('Name','Path A: Injector Characterisation','Position',[100 100 1100 620]);

% Plot 1: Time trace vs model predictions
subplot(2,2,1);
plot(t_trim, mdot*1000, 'k-', 'LineWidth',2, 'DisplayName','Measured N2O');
hold on;
for k = 1:n_d
    yline(mdot_SPI(k)*1000,  '--', 'Color',colors(k,:), 'LineWidth',1.5, ...
        'DisplayName',sprintf('SPI  D=%.4f"', inj.D_inj(k)));
    yline(mdot_Dyer(k)*1000, '-',  'Color',colors(k,:), 'LineWidth',2.0, ...
        'DisplayName',sprintf('Dyer D=%.4f"', inj.D_inj(k)));
end
yline(mdot_avg*1000, ':k', 'LineWidth',1.2, 'DisplayName',sprintf('Avg %.2f g/s', mdot_avg*1000));
xlabel('Time [s]'); ylabel('N2O Flow [g/s]');
title('Measured Flow vs Model Predictions','FontSize',11);
legend('Location','best','FontSize',7); grid on;

% Plot 2: Effective Cd vs design Cd per drill bit
subplot(2,2,2);
b = bar(1:n_d, [Cd_eff(:), repmat(Cd,n_d,1)], 0.6);
b(1).FaceColor = [0.2 0.5 0.8];
b(2).FaceColor = [0.82 0.82 0.82];
xticks(1:n_d);
xticklabels(arrayfun(@(d) sprintf('%.4f"',d), inj.D_inj, 'UniformOutput',false));
legend({'Effective Cd (from data)','Design Cd'}, 'Location','best','FontSize',9);
xlabel('Drill Bit Size'); ylabel('Discharge Coefficient');
title('Effective C_d vs Design C_d','FontSize',11); grid on;

% Plot 3: N2O flow vs downstream pressure sweep
subplot(2,2,3);
for k = 1:n_d
    plot(P_sweep_psi, mdot_Dyer_sw(k,:)*1000, '-',  'Color',colors(k,:), ...
        'LineWidth',2,   'DisplayName',sprintf('Dyer D=%.4f"',inj.D_inj(k)));
    hold on;
    plot(P_sweep_psi, mdot_SPI_sw(k,:)*1000,  '--', 'Color',colors(k,:), ...
        'LineWidth',1.2, 'HandleVisibility','off');
end
yline(mdot_avg*1000, '-.k', 'LineWidth',1.5, ...
    'DisplayName',sprintf('Measured avg %.2f g/s', mdot_avg*1000));
xline(14.7, ':k',  'LineWidth',1.0, 'DisplayName','Cold flow (atm)');
xline(440,  '--g', 'LineWidth',1.5, 'DisplayName','Target Pc = 440 psi');
xlabel('Downstream Pressure [psi]'); ylabel('N2O Flow [g/s]');
title('N2O Flow vs Downstream Pressure','FontSize',11);
legend('Location','northeast','FontSize',8); grid on;

% Plot 4: Two-phase correction magnitude (Dyer vs SPI)
subplot(2,2,4);
for k = 1:n_d
    pct = (mdot_Dyer_sw(k,:) - mdot_SPI_sw(k,:)) ./ mdot_Dyer_sw(k,:) * 100;
    plot(P_sweep_psi, pct, '-', 'Color',colors(k,:), 'LineWidth',2, ...
        'DisplayName',sprintf('D=%.4f"',inj.D_inj(k)));
    hold on;
end
yline(0,  '--k', 'LineWidth',1.0);
yline(10, ':r',  'LineWidth',1.5, 'DisplayName','10% significance threshold');
xline(14.7, ':k',  'LineWidth',1.0, 'HandleVisibility','off');
xline(440,  '--g', 'LineWidth',1.5, 'HandleVisibility','off');
xlabel('Downstream Pressure [psi]'); ylabel('Dyer - SPI [%]');
title('Two-Phase Correction Magnitude','FontSize',11);
legend('Location','best','FontSize',8); grid on;

sgtitle('Path A  --  N2O Injector Characterisation','FontSize',14,'FontWeight','bold');
end

% ---------------------------------------------------------
function targets = get_target_conditions()
% Collect N2O system conditions for injector sizing (Path B).

fprintf('\n── SYSTEM CONDITIONS ───────────────────────────\n');
fprintf('  Press ENTER to accept pressure budget defaults.\n\n');

targets.P_inlet      = parse_default(strtrim(input('  P_inlet      [psi]  (default 577):  ','s')), 577);
targets.P_downstream = parse_default(strtrim(input('  P_downstream [psi]  (default 14.7): ','s')), 14.7);
targets.P_hotfire    = parse_default(strtrim(input('  Target Pc    [psi]  (default 440):  ','s')), 440);
targets.Cd           = parse_default(strtrim(input('  Cd                  (default 0.7):  ','s')), 0.7);

fprintf('\n  P_inlet=%.1f psi  |  P_downstream=%.1f psi\n', ...
    targets.P_inlet, targets.P_downstream);
fprintf('  Target Pc=%.0f psi  |  Cd=%.2f\n', targets.P_hotfire, targets.Cd);
end

% ---------------------------------------------------------
function injector_sizing(mdot_avg, mdot_peak, targets)
% PATH B: Recommend drill bit sizes from measured N2O flow.
%
% Two conditions are evaluated for every diameter:
%   CF = cold flow  (downstream = atmospheric, no chamber)
%   HF = hot fire   (downstream = target chamber pressure)
%
% The sizing table flags each drill bit as:
%   * MATCH  -- within 10% of your measured avg flow at hot fire
%   ~ OK     -- reasonable but not optimal
%   ^ HIGH   -- would flow more than your measured peak
%   v LOW    -- would flow less than 85% of your measured avg

fprintf('\n── PATH B: INJECTOR SIZING ─────────────────────\n');

rho_l = 793;  rho_g = 128;  h_fg = 323.6e3;  P_sat = 5169e3;

P_inlet_pa = targets.P_inlet      * 6894.76;
P_cf_pa    = targets.P_downstream * 6894.76;
P_hf_pa    = targets.P_hotfire    * 6894.76;
Cd         = targets.Cd;

    function mA = dyer_mpa(dP)
        % mdot per unit area using Dyer two-phase model [kg/s/m2]
        if dP <= 0, mA = 0; return; end
        spi     = Cd * sqrt(2 * rho_l * dP);
        x       = min(max(dP/(rho_l*h_fg), 0), 1);
        rho_2ph = 1/(x/rho_g + (1-x)/rho_l);
        hem     = Cd * rho_2ph * sqrt(2*dP/rho_l);
        if P_inlet_pa >= P_sat, k = 0;
        else, k = sqrt(min(max(dP/max(P_inlet_pa-P_sat,1),0),1)); end
        mA = (1-k)*spi + k*hem;
    end

mpa_cf     = dyer_mpa(P_inlet_pa - P_cf_pa);
mpa_hf     = dyer_mpa(P_inlet_pa - P_hf_pa);
mpa_spi_cf = Cd * sqrt(2*rho_l*max(P_inlet_pa-P_cf_pa,0));
mpa_spi_hf = Cd * sqrt(2*rho_l*max(P_inlet_pa-P_hf_pa,0));

% Back-calculate diameter from measured mdot
fprintf('\n  ── Diameter matching your measured flow ──\n');
labels = {'Average mdot','Peak mdot'};
mvals  = [mdot_avg, mdot_peak];
for i = 1:2
    D_cf = sqrt(4*(mvals(i)/mpa_cf)/pi) / 0.0254;
    D_hf = sqrt(4*(mvals(i)/mpa_hf)/pi) / 0.0254;
    fprintf('  %s (%.4f kg/s):\n', labels{i}, mvals(i));
    fprintf('    Cold flow downstream (%.1f psi):  D = %.4f in\n', targets.P_downstream, D_cf);
    fprintf('    Hot fire Pc (%.0f psi):           D = %.4f in\n\n', targets.P_hotfire, D_hf);
end

% Continuous diameter sweep for plots
D_sweep     = linspace(0.04, 0.22, 300);
A_sweep     = pi/4 .* (D_sweep*0.0254).^2;
mdot_cf     = mpa_cf    .* A_sweep;
mdot_hf     = mpa_hf    .* A_sweep;
mdot_spi_cf = mpa_spi_cf .* A_sweep;
mdot_spi_hf = mpa_spi_hf .* A_sweep;

% Best-fit diameter at hot fire conditions
A_best  = mdot_avg / mpa_hf;
D_best  = sqrt(4*A_best/pi) / 0.0254;

% Flow vs downstream pressure for best-fit diameter
P_sweep_pa  = linspace(P_cf_pa, P_hf_pa, 300);
P_sweep_psi = P_sweep_pa / 6894.76;
mdot_Pbest  = arrayfun(@(p) dyer_mpa(P_inlet_pa-p)*A_best, P_sweep_pa);

% Sizing table: standard drill bit sizes
D_table     = [0.052, 0.063, 0.070, 0.076, 0.086, 0.094, 0.098, 0.104, 0.110, 0.120, 0.136, 0.147];
drill_names = {'#55','#52','#50','#48','#44','#42','#40','#38','#35','#31','#29','#26'};

fprintf('  ┌──────────┬───────────┬───────────┬───────────┬──────────────┐\n');
fprintf('  │ D [in]   │ mdot CF   │ mdot HF   │ Note      │ Drill bit    │\n');
fprintf('  │          │ [g/s]Dyer │ [g/s]Dyer │           │              │\n');
fprintf('  ├──────────┼───────────┼───────────┼───────────┼──────────────┤\n');
for ki = 1:length(D_table)
    d    = D_table(ki);
    A_d  = pi/4*(d*0.0254)^2;
    m_cf = mpa_cf * A_d * 1000;
    m_hf = mpa_hf * A_d * 1000;
    if     m_hf > mdot_peak*1000*1.05,                     note = '^ HIGH  ';
    elseif m_hf < mdot_avg*1000*0.85,                      note = 'v LOW   ';
    elseif abs(m_hf-mdot_avg*1000) < mdot_avg*1000*0.10,   note = '* MATCH ';
    else,                                                   note = '~ OK    '; end
    cur = '';
    if abs(d-0.094) < 0.001, cur = ' <- current'; end
    fprintf('  │ D=%.4f │ %7.3f   │ %7.3f   │ %s  │ %-6s%-6s│\n', ...
        d, m_cf, m_hf, note, drill_names{ki}, cur);
end
fprintf('  └──────────┴───────────┴───────────┴───────────┴──────────────┘\n');
fprintf('  Measured avg=%.2f g/s,  peak=%.2f g/s\n', mdot_avg*1000, mdot_peak*1000);
fprintf('  CF=cold flow (atm),  HF=hot fire Pc=%.0f psi\n\n', targets.P_hotfire);

figure('Name','Path B: Injector Sizing','Position',[100 100 1100 620]);

% Plot 1: N2O flow vs orifice diameter
subplot(2,2,1);
plot(D_sweep, mdot_cf*1000,     'b-',  'LineWidth',2.5, 'DisplayName','Cold flow (Dyer)');
hold on;
plot(D_sweep, mdot_hf*1000,     'r-',  'LineWidth',2.5, 'DisplayName','Hot fire (Dyer)');
plot(D_sweep, mdot_spi_cf*1000, 'b--', 'LineWidth',1.2, 'DisplayName','Cold flow (SPI)');
plot(D_sweep, mdot_spi_hf*1000, 'r--', 'LineWidth',1.2, 'DisplayName','Hot fire (SPI)');
yline(mdot_avg*1000,  '-.', 'Color',[1 0.5 0], 'LineWidth',1.8, ...
    'DisplayName',sprintf('Meas avg %.2f g/s',  mdot_avg*1000));
yline(mdot_peak*1000, ':',  'Color',[1 0.5 0], 'LineWidth',1.8, ...
    'DisplayName',sprintf('Meas peak %.2f g/s', mdot_peak*1000));
xline(D_best, '-.m', 'LineWidth',1.5, 'DisplayName',sprintf('Best-fit D=%.4f in', D_best));
xlabel('Orifice Diameter [in]'); ylabel('N2O Mass Flow [g/s]');
title('N2O Flow vs Orifice Diameter','FontSize',11);
legend('Location','northwest','FontSize',8); grid on;

% Plot 2: N2O flow vs downstream pressure (best-fit diameter)
subplot(2,2,2);
plot(P_sweep_psi, mdot_Pbest*1000, 'r-', 'LineWidth',2.5, ...
    'DisplayName',sprintf('D=%.4f in (best-fit)', D_best));
hold on;
yline(mdot_avg*1000,  '--b', 'LineWidth',1.8, ...
    'DisplayName',sprintf('Meas avg %.2f g/s',  mdot_avg*1000));
yline(mdot_peak*1000, ':b',  'LineWidth',1.8, ...
    'DisplayName',sprintf('Meas peak %.2f g/s', mdot_peak*1000));
xline(targets.P_downstream, ':k',  'LineWidth',1.2, ...
    'DisplayName',sprintf('Cold flow %.1f psi', targets.P_downstream));
xline(targets.P_hotfire,    '--g', 'LineWidth',1.5, ...
    'DisplayName',sprintf('Hot fire Pc=%.0f psi', targets.P_hotfire));
xlabel('Downstream Pressure [psi]'); ylabel('N2O Flow [g/s]');
title(sprintf('N2O Flow vs Downstream Pressure  (D=%.4f in)', D_best),'FontSize',11);
legend('Location','northeast','FontSize',8); grid on;

% Plot 3: SPI vs Dyer correction at hot fire conditions
subplot(2,2,3);
pct_diff = (mdot_hf - mdot_spi_hf) ./ mdot_hf * 100;
plot(D_sweep, pct_diff, 'm-', 'LineWidth',2);
hold on;
yline(0,  '--k', 'LineWidth',1.0);
yline(10, ':r',  'LineWidth',1.5, 'DisplayName','10% significance');
xline(D_best, '-.m', 'LineWidth',1.5, 'DisplayName',sprintf('Best-fit D=%.4f in', D_best));
xlabel('Orifice Diameter [in]'); ylabel('Dyer - SPI [%]');
title('Two-Phase Correction: Dyer vs SPI  (Hot Fire cond.)','FontSize',11);
legend('Location','best','FontSize',9); grid on;

% Plot 4: Injector pressure ratio stability check
subplot(2,2,4);
PR_val = targets.P_inlet / targets.P_hotfire;
fill([min(D_sweep), max(D_sweep), max(D_sweep), min(D_sweep)], ...
     [0, 0, 1.2, 1.2], [1 0.85 0.85], 'EdgeColor','none', ...
     'DisplayName','Unstable zone (PR < 1.2)');
hold on;
plot(D_sweep, repmat(PR_val, size(D_sweep)), 'r-', 'LineWidth',2.5, ...
    'DisplayName',sprintf('System PR = %.3f  (P_{inlet}/P_c)', PR_val));
yline(1.2, '--r', 'LineWidth',1.5, 'DisplayName','Min stable PR = 1.2');
xline(D_best, '-.m', 'LineWidth',1.5, ...
    'DisplayName',sprintf('Best-fit D=%.4f in', D_best));
xlabel('Orifice Diameter [in]'); ylabel('Pressure Ratio  P_{inlet}/P_c');
title('Injector Pressure Ratio  (Stability Check)','FontSize',11);
legend('Location','best','FontSize',8); grid on;
ylim([0.8, PR_val*1.4]);

sgtitle('Path B  --  N2O Injector Sizing','FontSize',14,'FontWeight','bold');

fprintf('  RECOMMENDATION:\n');
fprintf('  Best-fit orifice:  D = %.4f in\n', D_best);
fprintf('  (Dyer model, matches avg N2O flow at Pc=%.0f psi, Cd=%.2f)\n', ...
    targets.P_hotfire, Cd);
fprintf('  System PR = %.3f  ', PR_val);
if PR_val >= 1.2
    fprintf('--  OK (above 1.2 stability threshold)\n');
else
    fprintf('--  WARNING: below 1.2, flow may be pressure-driven by chamber\n');
end
fprintf('  If Dyer-SPI correction > 10%%: use Dyer model for final sizing.\n\n');
end

% ---------------------------------------------------------
function val = parse_default(str, default)
% Return default if input string is empty, else parse as number.
if isempty(str), val = default;
else,            val = str2double(str); end
end