function CB = constant()
% Centralized constants & defaults for HALCS-Net

%% Moon and Earth geometry and Constants

CB.moon.name        = "Moon";
CB.moon.R_km        = 1737.4;              % mean radius

CB.moon.R_m         = CB.moon.R_km * 1e3;
CB.eath_moon.R      = 384400;

CB.k_dB = -228.6;
CB.C = 299792458;

%% L2 geometry

CB.l2.R_L2_km       = 61350; % distance of moon to L2
CB.l2.Ay_km         = 39000;
CB.l2.Az_km         = 20000;
CB.links.RL.range_km = [63000 71000];

%% Frequencies allocation

CB.ka.freq_GHz.RL_UL = 27.25;
CB.ka.freq_GHz.RL_DL = 23.35;
CB.ka.freq_GHz.RE_UL = 23.00;
CB.ka.freq_GHz.RE_DL = 26.25;
CB.ka.freq_GHz.LE_UL = 27.25;
CB.ka.freq_GHz.LE_DL = 26.25;

%% Relay, UEs and Ground stations antenna's properties

% Antenna efficiencies
CB.ant.eta   = 0.60;

% Relay antenna property
CB.relay.moonBeam.HPBW_deg  = 3.07;
CB.relay.earthBeam.HPBW_deg = 1.632;

% Relay powers
CB.relay.moon.Ptx_ref_W = 40;
CB.relay.earth.Ptx_ref_W = 100;

% UE anttena
CB.ue.dish_ref_m = 0.7;
CB.ue.dish_m     = CB.ue.dish_ref_m - 0.2 : 0.05 : CB.ue.dish_ref_m + 0.2; % simulate 20cm below and above 0.7
CB.ue.Ptx_ref_W  = 20; % power of UE
CB.ue.Ptx_W      = [20 30 40]; % power range for simulation

CB.ue.earthBeam.HPBW_deg = 1.11;

% Ground station antenna's properties
CB.gs.dish_m     = [8.2 11 18 34];
CB.gs.dish_ref     = 34;
CB.gs.Ptx_W      = 500;

%% Loss and Noise

% Minimal loss
CB.loss.pol_dB    = 0.8; % we can ignore this
CB.loss.point_dB  = 0.5;
CB.loss.miscTx_dB = CB.loss.point_dB + CB.loss.pol_dB;

CB.loss.groundK = 20;

% Noise temperatures
CB.noise.rxK_relay_RL = 401; 
CB.noise.rxK_relay_RE = 400;

CB.noise.rxK_ue_RL    = [110 170];
CB.noise.rxK_ue_RE    = 400;

CB.noise.gsNF_dB      = 1.0;

%% Links properties

CB.link.Rb_bps = 10e6; % 10mbps
CB.link.min_EbN0 = 2; % coded QPSK

%% Ground sites

CB.sites.catalog = struct( ...
  'name', {'Chajnantor','Teide','Hanle','Ngari', 'Goldstone','INO' ...
            }, ...
  'lla',  { [-23.02,  -67.76, 5.05], ...   % Chajnantor
            [ 28.27,  -16.64, 2.39], ...   % Teide
            [ 32.78,   78.96, 4.50], ...   % Hanle
            [ 32.32,   80.03, 5.10], ...   % Ngari
            [ 35.247,-116.793, 1.01], ...  % Goldstone
            [ 33.674,  51.318, 3.57], ...  % INO
          });

%% Simulation defaults

CB.defaults.minElevation_deg = 20; % minimum elevation
end