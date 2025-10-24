%% L2 simulation of UE <-> Relay

clear; clc;
CB = constant();

%% Helper functions and numbers
c                     = CB.C;
lambda                = @(f_GHz) c/(f_GHz*1e9); % getting wave lenght
dishGain_dBi          = @(D,eta,f_GHz) 10*log10(eta*(pi*D/lambda(f_GHz)).^2); %getting dishgain
gainFromHPBW_dBi      = @(HPBW_deg,eta) 10*log10(eta*41253/(HPBW_deg.^2)); %getting dishgain from HPBW
FSPL_dB               = @(rng_km,f_GHz) 92.45 + 20*log10(f_GHz) + 20*log10(rng_km); %calculation of Free space loss
W2dBW                 = @(P_W) 10*log10(P_W); % watt to dbW
k_dB                  = CB.k_dB;
Rb_dBHz               = 10*log10(CB.link.Rb_bps); % data rate to dB

% Relay gain
Grelay_moon_dBi       = gainFromHPBW_dBi(CB.relay.moonBeam.HPBW_deg, CB.ant.eta);

% Frequencies, ranges
fUL                   = CB.ka.freq_GHz.RL_UL;
fDL                   = CB.ka.freq_GHz.RL_DL;
rng_near_km           = CB.links.RL.range_km(1);
rng_far_km            = CB.links.RL.range_km(2);

% Efficiencies & losses
eta                = CB.ant.eta;
Ltx_dB                = CB.loss.miscTx_dB;

% UE and Relay Antenna's property
D_UE_ref_m            = CB.ue.dish_ref_m;
P_UE_ref_W            = CB.ue.Ptx_ref_W;
P_Relay_ref_W         = CB.relay.moon.Ptx_ref_W;

% Numbers for range simulation
D_UE_list_m           = CB.ue.dish_m;
P_UE_list_W           = CB.ue.Ptx_W;
P_Relay_list_W        = [20 30 40];

% Temperatures
Trelay_lo_hi_K        = CB.noise.rxK_relay_RL;
Tue_lo_hi_K           = [CB.noise.rxK_ue_RL(1),    CB.noise.rxK_ue_RL(2)];
Trelay_hi_K           = CB.noise.rxK_relay_RE;
Tue_hi_K              = CB.noise.rxK_ue_RL(2);

%% Calculations

% UL (from UE to Relay)
ul_ebn0 = @(P_UE_W, D_UE_m, TrelayK, rng_km) ...
    ( W2dBW(P_UE_W) + dishGain_dBi(D_UE_m, eta, fUL) ...
      - Ltx_dB ...
      + (Grelay_moon_dBi - 10*log10(TrelayK)) ...
      - FSPL_dB(rng_km, fUL) - k_dB ) - Rb_dBHz;

% DL (from Relay to UE)
dl_ebn0 = @(P_R_W, D_UE_m, TueK, rng_km) ...
    ( W2dBW(P_R_W) + Grelay_moon_dBi ...
      - Ltx_dB ...
      + (dishGain_dBi(D_UE_m, eta, fDL) - 10*log10(TueK)) ...
      - FSPL_dB(rng_km, fDL) - k_dB ) - Rb_dBHz;

%% Eb/N0 vs UE dish at FAR range, HIGHEST temperature

cols = lines(numel(P_UE_list_W));
Eb_UL_D = zeros(numel(P_UE_list_W), numel(D_UE_list_m));
for i=1:numel(P_UE_list_W)
    Eb_UL_D(i,:) = arrayfun(@(D) ul_ebn0(P_UE_list_W(i), D, Trelay_hi_K, rng_far_km), D_UE_list_m);
end

Eb_DL_D = zeros(numel(P_Relay_list_W), numel(D_UE_list_m));
for j=1:numel(P_Relay_list_W)
    Eb_DL_D(j,:) = arrayfun(@(D) dl_ebn0(P_Relay_list_W(j), D, Tue_hi_K, rng_far_km), D_UE_list_m);
end

% scale from min -1 to max+2
ymax2 = max([Eb_UL_D(:); Eb_DL_D(:)]);
ymax2 = ceil(ymax2*2)/2 + 2;
ymin2 = min([Eb_UL_D(:); Eb_DL_D(:)]) - 1;

xmax2 = max(D_UE_list_m(:)) + 0.01;
xmin = min(D_UE_list_m(:)) - 0.01;

figure('Color','w','Name','RL: Eb/N0 vs UE dish (farthest range, worst Tsys)');
hold on; grid on; box on;

% UL
for i=1:numel(P_UE_list_W)
    plot(D_UE_list_m, Eb_UL_D(i,:), 'o-','LineWidth',1.8,'Color',cols(i,:), ...
         'DisplayName', sprintf('UL  P_{UE}=%d W', P_UE_list_W(i)));
end

% DL (reuse colors but with dash)
for j=1:numel(P_Relay_list_W)
    plot(D_UE_list_m, Eb_DL_D(j,:), 's--','LineWidth',1.8,'Color',cols(j,:), ...
         'DisplayName', sprintf('DL  P_{Relay}=%d W', P_Relay_list_W(j)));
end

yline(CB.link.min_EbN0,':', string(CB.link.min_EbN0) + ' dB (coded QPSK)','LabelHorizontalAlignment','left');
xlabel('UE dish diameter (m)'); ylabel('E_b/N_0 (dB)');
title(sprintf('E_b/N_0 vs UE dish — RL link  (far=%.0f km, worst T_{relay}/T_{UE})', rng_far_km));
legend('Location','best');
ylim([ymin2 ymax2]);
xlim([xmin xmax2]);