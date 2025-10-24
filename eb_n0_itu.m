%% Simulaiton of Earth <-> UE and Relay

clear; clc;
CB = constant();

%% Ensure ITU digital maps for P.618 are available
have_maps = exist('maps.mat','file') && exist('p836.mat','file') && ...
            exist('p837.mat','file') && exist('p840.mat','file');
if ~have_maps
    try
        websave('ITURDigitalMaps.tar.gz', ...
            'https://www.mathworks.com/supportfiles/spc/P618/ITURDigitalMaps.tar.gz');
        untar('ITURDigitalMaps.tar.gz');
        addpath(pwd);
    catch ME
        warning('Could not download ITU digital maps: %s.', ME.message);
    end
end

%% Helper functions and numbers
c                     = CB.C;
lambda                = @(f_GHz) c/(f_GHz*1e9); % getting wave lenght
dishGain_dBi          = @(D,eta,f_GHz) 10*log10(eta*(pi*D/lambda(f_GHz)).^2); %getting dishgain
gainFromHPBW_dBi      = @(HPBW_deg,eta) 10*log10(eta*41253/(HPBW_deg.^2)); %getting dishgain from HPBW
FSPL                  = @(rng_km,f_GHz) 92.45 + 20*log10(f_GHz) + 20*log10(rng_km); %calculation of Free space loss
W2dBW                 = @(P_W) 10*log10(P_W); % watt to dbW
k_dB                  = CB.k_dB;
Rb_dBHz               = 10*log10(CB.link.Rb_bps); % data rate to dB

% Relay gain
G_relay_earth_dBi = gainFromHPBW_dBi(CB.relay.earthBeam.HPBW_deg, CB.ant.eta);

% TX losses
Ltx_misc_dB = CB.loss.miscTx_dB;

% Sites
sitesCat = CB.sites.catalog;
sites = sitesCat;
Ns = numel(sites);

% list of dishes
gsD_list = CB.gs.dish_m;

% Minimum elevation of ground stations
elev_deg = 15;

D_EM   = CB.eath_moon.R;
D_EL2  = D_EM + CB.l2.R_L2_km; %L2 from earth
r_perp = hypot(CB.l2.Ay_km, CB.l2.Az_km); % max distance possible of l2 satellite
far_RE_km = sqrt(D_EL2^2 + r_perp^2); % max distance of l2 satellite to earth
far_LE_km = D_EM; % UE -> Earth distance

% Availability choices
pGrid = [1 0.05 0.01 0.005 0.001];

% Threshold is the same as the qpsk coded but for margin
EbN0_thr_dB = CB.link.min_EbN0;

clrUL = [0.20 0.60 1.00];
clrDL = [1.00 0.55 0.10]; 
lsUL  = '-';  lsDL = '--';

%% RE (Relay <-> Earth)

fRE = figure('Color','w','Name','RE: Eb/N0 vs Availability');
tlRE  = tiledlayout(fRE,2,2,'TileSpacing','compact','Padding','compact');
title(tlRE, sprintf('Relay<->Earth, Elev=%d°, Range=%.1f Mm', ...
      elev_deg, far_RE_km/1000), 'FontWeight','bold');

hUL = gobjects(Ns, 1);
hDL = gobjects(Ns, 1);

colors = lines(Ns);

for ig = 1:numel(gsD_list)
    Dgs = gsD_list(ig);

    % EIRPs for atmoshpere
    EIRP_UL = W2dBW(CB.gs.Ptx_W)              + dishGain_dBi(Dgs, CB.ant.eta, CB.ka.freq_GHz.RE_UL)  - Ltx_misc_dB;
    EIRP_DL = W2dBW(CB.relay.earth.Ptx_ref_W) + G_relay_earth_dBi                                   - Ltx_misc_dB;

    % Relay RX side GT
    GT_relay_RX = G_relay_earth_dBi - 10*log10(CB.noise.rxK_relay_RE);

    nexttile; hold on; grid on; box on;
    set(gca,'XScale','log'); xlim([min(pGrid) max(pGrid)]); set(gca,'XDir','reverse');
    xlabel('Availability exceedance p (%)'); ylabel('E_b/N_0 (dB)');
    title(sprintf('GS dish = %.1f m', Dgs));

    % Per-site calculation
    for s = 1:Ns
        LLA = sites(s).lla;

        Eb_UL = nan(size(pGrid));
        Eb_DL = nan(size(pGrid));

        for ip = 1:numel(pGrid)
            av = struct('name',sprintf('p=%.5g%%',pGrid(ip)), ... % for propagation function
                        'TotalAnnualExceedance', pGrid(ip), ...
                        'RainAnnualExceedance',  pGrid(ip));

            % P.618 calculation
            [GT_gs_DL, L_atm_DL] = gs_GT_and_L(LLA, Dgs, CB.ant.eta, elev_deg, CB.ka.freq_GHz.RE_DL, CB.noise.gsNF_dB, av, CB.loss.groundK);
            [~,        L_atm_UL] = gs_GT_and_L(LLA, Dgs, CB.ant.eta, elev_deg, CB.ka.freq_GHz.RE_UL, CB.noise.gsNF_dB, av, 0);

            % UL Ground staion -> Relay
            Eb_UL(ip) = EIRP_UL + GT_relay_RX - FSPL(far_RE_km, CB.ka.freq_GHz.RE_UL) - (k_dB + L_atm_UL) - Rb_dBHz;
            % DL Relay to Ground stations
            Eb_DL(ip) = EIRP_DL + GT_gs_DL    - FSPL(far_RE_km, CB.ka.freq_GHz.RE_DL) - (k_dB + L_atm_DL) - Rb_dBHz;
        end

        hUL(s) = plot(pGrid, Eb_UL, lsUL, 'LineWidth', 1.6, 'Color', colors(s,:));
        hDL(s) = plot(pGrid, Eb_DL, lsDL, 'LineWidth', 1.6, 'Color', colors(s,:));
    end

    yline(EbN0_thr_dB,':', sprintf('%.1f dB (coded QPSK)',EbN0_thr_dB), 'LabelHorizontalAlignment','left');
end

% Build labels on the center of sub plots
axForLegend = gca;                         % get last tile
names  = string({sites.name}).';
labels = [names + " UL"; names + " DL"];
allH   = [hUL(:); hDL(:)];

nItems = numel(allH);
nCols  = 6;

lgd = legend(axForLegend, allH, labels, ...
    'Location','southoutside', ...
    'NumColumns', nCols, ...
    'Interpreter','none', ...
    'FontSize', 9, ...
    'Box','off');

lgd.Layout.Tile = 'south'; 


%% LE (UE <-> Earth)

D_UE      = CB.ue.dish_ref_m;
T_UE_K    = CB.noise.rxK_ue_RE;
GT_UE_RX  = dishGain_dBi(D_UE, CB.ant.eta, CB.ka.freq_GHz.LE_DL) - 10*log10(T_UE_K);

fLE = figure('Color','w','Name','LE: Eb/N0 vs Availability');
tlLE  = tiledlayout(fLE,2,2,'TileSpacing','compact','Padding','compact');
title(tlLE, sprintf('UE<->Earth, Elev=%d° (GS), Range=%.1f Mm', ...
      elev_deg, far_LE_km/1000), 'FontWeight','bold');

for ig = 1:numel(gsD_list)
    Dgs = gsD_list(ig);

    % EIRPs
    EIRP_Down = W2dBW(CB.ue.Ptx_ref_W) + dishGain_dBi(D_UE, CB.ant.eta, CB.ka.freq_GHz.LE_UL) - Ltx_misc_dB;
    EIRP_Up   = W2dBW(CB.gs.Ptx_W)     + dishGain_dBi(Dgs,  CB.ant.eta, CB.ka.freq_GHz.LE_DL) - Ltx_misc_dB;

    nexttile; hold on; grid on; box on;
    set(gca,'XScale','log'); xlim([min(pGrid) max(pGrid)]); set(gca,'XDir','reverse');
    xlabel('Availability exceedance p (%)'); ylabel('E_b/N_0 (dB)');
    title(sprintf('GS dish = %.1f m  (UE dish = %.2f m)', Dgs, D_UE));

    for s = 1:Ns
        LLA = sites(s).lla;

        Eb_Down = nan(size(pGrid));
        Eb_Up   = nan(size(pGrid));

        for ip = 1:numel(pGrid)
            av = struct('name',sprintf('p=%.5g%%',pGrid(ip)), ... % for atmoshpher funciton
                        'TotalAnnualExceedance', pGrid(ip), ...
                        'RainAnnualExceedance',  pGrid(ip));

            % P.618 calculation
            [GT_gs_Down, L_atm_Down] = gs_GT_and_L(LLA, Dgs, CB.ant.eta, elev_deg, CB.ka.freq_GHz.LE_UL, CB.noise.gsNF_dB, av, CB.loss.groundK);
            [~,          L_atm_Up]   = gs_GT_and_L(LLA, Dgs, CB.ant.eta, elev_deg, CB.ka.freq_GHz.LE_DL, CB.noise.gsNF_dB, av, 0);

            % UE -> Earth (downlink)
            Eb_Down(ip) = EIRP_Down + GT_gs_Down - FSPL(far_LE_km, CB.ka.freq_GHz.LE_UL) - (k_dB + L_atm_Down) - Rb_dBHz;

            % Earth -> UE (uplink)
            Eb_Up(ip)   = EIRP_Up   + GT_UE_RX  - FSPL(far_LE_km, CB.ka.freq_GHz.LE_DL) - (k_dB + L_atm_Up)   - Rb_dBHz;
        end

        hUL(s) = plot(pGrid, Eb_Up,   lsUL, 'LineWidth',1.6, 'Color', colors(s,:)); % Uplink (from Earth)
        hDL(s) = plot(pGrid, Eb_Down, lsDL, 'LineWidth',1.6, 'Color', colors(s,:)); % Downlink (to Earth)
    end

    yline(EbN0_thr_dB,':', sprintf('%.1f dB (coded QPSK)',EbN0_thr_dB), 'LabelHorizontalAlignment','left');
end

axForLegend = gca;                         % get last tile
names  = string({sites.name}).';
labels = [names + " UL"; names + " DL"];
allH   = [hUL(:); hDL(:)];

nItems = numel(allH);
nCols  = 6;

lgd = legend(axForLegend, allH, labels, ...
    'Location','southoutside', ...
    'NumColumns', nCols, ...
    'Interpreter','none', ...
    'FontSize', 9, ...
    'Box','off');

lgd.Layout.Tile = 'south'; 

%% Margin plot

mRE_UL = zeros(1,Ns);  mRE_DL = zeros(1,Ns);
Dgs_ref = CB.gs.dish_ref;
pNeedle = 0.001;

for s=1:Ns
    LLA = sites(s).lla;
    [mRE_UL(s), mRE_DL(s)] = compute_margin(Dgs_ref, LLA, CB, ...
                                  CB.ka.freq_GHz.RE_UL, CB.ka.freq_GHz.RE_DL, CB.noise.rxK_relay_RE, ...
                                  CB.relay.earthBeam.HPBW_deg, CB.relay.earth.Ptx_ref_W, ...
                                  far_RE_km, EbN0_thr_dB, elev_deg, pNeedle);
end

plot_margins({sites.name}, mRE_UL, mRE_DL, ...
    sprintf('Margins @ p=0.001%% (Relay<->Earth, GS=%.0f m, Range=%.1f Mm)', Dgs_ref, far_RE_km/1000));

mLE_Up   = zeros(1,Ns);
mLE_Down = zeros(1,Ns);
for s=1:Ns
    LLA = sites(s).lla;
    [mLE_Up(s), mLE_Down(s)] = compute_margin(Dgs_ref, LLA, CB, ...
                                  CB.ka.freq_GHz.LE_UL, CB.ka.freq_GHz.LE_DL, CB.noise.rxK_ue_RE, ...
                                  CB.ue.earthBeam.HPBW_deg, CB.ue.Ptx_ref_W, ...
                                  far_LE_km, EbN0_thr_dB, elev_deg, pNeedle);
end

plot_margins({sites.name}, mLE_Up, mLE_Down, ...
    sprintf('Margins @ p=0.001%% (UE<->Earth, GS=%.0f m, Range=%.1f Mm)', Dgs_ref, far_LE_km/1000));

%% Helper functions
function [GT_dB, L_atm_dB] = gs_GT_and_L(siteLLA, D_rx_m, eta, elev_deg, f_GHz, gsNF_dB, avail, add_tground)
    cfgP = p618Config;
    cfgP.Frequency             = f_GHz*1e9;
    cfgP.ElevationAngle        = elev_deg;
    cfgP.Latitude              = siteLLA(1);
    cfgP.Longitude             = siteLLA(2);
    cfgP.AntennaDiameter       = D_rx_m;
    cfgP.AntennaEfficiency     = eta;
    cfgP.TotalAnnualExceedance = avail.TotalAnnualExceedance;
    cfgP.RainAnnualExceedance  = avail.RainAnnualExceedance;

    [pl, ~, T_sky] = p618PropagationLosses(cfgP);
    L_atm_dB = pl.At;

    G_rec_dBi = 10*log10( eta*(pi*D_rx_m/(299792458/(f_GHz*1e9)))^2 );
    T_rxK    = 290*(10^(gsNF_dB/10)-1);
    T_sysK   = T_sky + T_rxK + add_tground;
    GT_dB = G_rec_dBi - 10*log10(T_sysK);
end

% computing marging for RE
function [mUL, mDL] = compute_margin(Dgs, LLA, CB, UP_Freq, Down_Freq, rxK, HPBW_deg, Pw, range_km, EbN0_thr_dB, elev_deg, pPercent)
    c                     = CB.C;
    lambda                = @(f_GHz) c/(f_GHz*1e9); % getting wave lenght
    dishGain_dBi          = @(D,eta,f_GHz) 10*log10(eta*(pi*D/lambda(f_GHz)).^2); %getting dishgain
    gainFromHPBW_dBi      = @(HPBW_deg,eta) 10*log10(eta*41253/(HPBW_deg.^2)); %getting dishgain from HPBW
    FSPL                  = @(rng_km,f_GHz) 92.45 + 20*log10(f_GHz) + 20*log10(rng_km); %calculation of Free space loss
    W2dBW                 = @(P_W) 10*log10(P_W); % watt to dbW
    k_dB                  = CB.k_dB;
    Rb_dBHz               = 10*log10(CB.link.Rb_bps); % data rate to dB
    Ltx = CB.loss.miscTx_dB;

    av = struct('TotalAnnualExceedance',pPercent, 'RainAnnualExceedance',pPercent);

    [~, L_atm_UL]  = gs_GT_and_L(LLA, Dgs, CB.ant.eta, elev_deg, UP_Freq, CB.noise.gsNF_dB, av, 0);
    EIRP           = W2dBW(CB.gs.Ptx_W) + dishGain_dBi(Dgs, CB.ant.eta, UP_Freq) - Ltx;
    Gt             = gainFromHPBW_dBi(HPBW_deg, CB.ant.eta);
    GT_RX          = Gt - 10*log10(rxK);
    Eb_UL          = EIRP + GT_RX - FSPL(range_km, UP_Freq) - (k_dB + L_atm_UL) - Rb_dBHz;
    mUL            = Eb_UL - EbN0_thr_dB;

    [GT_gs_DL, L_atm_DL] = gs_GT_and_L(LLA, Dgs, CB.ant.eta, elev_deg, Down_Freq, CB.noise.gsNF_dB, av, CB.loss.groundK);
    EIRP                 = W2dBW(Pw) + gainFromHPBW_dBi(HPBW_deg, CB.ant.eta) - Ltx;
    Eb_DL                = EIRP + GT_gs_DL - FSPL(range_km, Down_Freq) - (k_dB + L_atm_DL) - Rb_dBHz;
    mDL                  = Eb_DL - EbN0_thr_dB;
end

function plot_margins(siteNames, mUL, mDL, ttl)
    labUL='Uplink (from Earth)';
    labDL='Downlink (to Earth)';
    N = numel(siteNames); x = 1:N;

    clrUL = [0.20 0.60 1.00];  % blue
    clrDL = [1.00 0.55 0.10];  % orange

    figure('Color','w','Name',ttl); hold on; grid on; box on;
    for i=1:N
        plot([x(i) x(i)], [0 mUL(i)], 'Color',clrUL,'LineWidth',2);
        plot(x(i), mUL(i), 'o', 'MarkerFaceColor',clrUL, 'MarkerEdgeColor','k', 'MarkerSize',6);
        plot([x(i)+0.25 x(i)+0.25], [0 mDL(i)], 'Color',clrDL,'LineWidth',2);
        plot(x(i)+0.25, mDL(i), 'o', 'MarkerFaceColor',clrDL, 'MarkerEdgeColor','k', 'MarkerSize',6);
    end
    xlim([0.5 N+0.75]); ylabel('Margin: E_b/N_0 - threshold (dB)');
    set(gca,'XTick',x+0.125,'XTickLabel',siteNames,'XTickLabelRotation',35);

    hUL = plot(nan, nan, '-o', 'Color', clrUL, 'LineWidth', 2, ...
               'MarkerFaceColor', clrUL, 'MarkerEdgeColor','k');
    hDL = plot(nan, nan, '-o', 'Color', clrDL, 'LineWidth', 2, ...
               'MarkerFaceColor', clrDL, 'MarkerEdgeColor','k');
    legend([hUL hDL], {labUL, labDL}, 'Location', 'northwest');
    title(ttl);
end
