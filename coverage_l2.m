clear; clc;

CB = constant();
Rm_km = CB.moon.R_km;
gridStep      = 2;
eMinEarth_deg = CB.defaults.minElevation_deg;
eMinRelay_deg = CB.defaults.minElevation_deg;

R_L2_km = CB.l2.R_L2_km;

Ay_km = CB.l2.Ay_km;
Az_km = CB.l2.Az_km;

r_perp   = Az_km;
alpha_deg = atan2d(r_perp, R_L2_km);
d_mean_km = hypot(R_L2_km, r_perp);

Nsat          = 4;
dL2_km        = [d_mean_km, d_mean_km, d_mean_km, d_mean_km];
phiOffsets_deg = 0:1:90;

[lonG, latG] = meshgrid(0:gridStep:360-gridStep, -90:gridStep:90);
clat = cosd(latG); slat = sind(latG);
clon = cosd(lonG); slon = sind(lonG);

% unit vec of each coordination
rx = Rm_km*clat.*clon;  ry = Rm_km*clat.*slon;  rz = Rm_km*slat;
nx = rx/Rm_km;          ny = ry/Rm_km;          nz = rz/Rm_km;

eEarth_deg = asind( max(min(nx,1),-1) );
EarthOK    = (eEarth_deg >= eMinEarth_deg);

pctAny    = zeros(size(phiOffsets_deg));
pctNone   = zeros(size(phiOffsets_deg));
worstR_km = zeros(size(phiOffsets_deg));
minest_km = zeros(size(phiOffsets_deg));
catMaps   = cell(size(phiOffsets_deg));

for ii = 1:numel(phiOffsets_deg)
    phi0 = phiOffsets_deg(ii);

    % building ring and compute the d vec
    [u_rel, dvec] = buildRing(Nsat, phi0, alpha_deg, dL2_km);
    Rpos = u_rel .* dvec;

    [res, ~, minRhoVis] = coverage_and_worst(Rpos, rx, ry, rz, nx, ny, nz, ...
                                             eMinRelay_deg, EarthOK, latG);
    pctAny(ii)  = res.pct_Any;
    pctNone(ii) = res.pct_None;
    worstR_km(ii) = maxfinite(minRhoVis(:));
    minest_km(ii) = minfinite(minRhoVis(:));
    catMaps{ii} = res.catMap;
end

[~, idxWorst]     = min(pctAny);
[worstR_all, iR]  = max(worstR_km);
[minR_all, iR]  = min(minest_km);

fprintf('Worst coverage %.2f %%  at offset=%d°\n', ...
        pctAny(idxWorst), phiOffsets_deg(idxWorst));
fprintf('Slant range and Shortest slant range: %.0f km  %.0f km\n\n', ...
        worstR_all, minR_all);

plotCategoryMap(catMaps{idxWorst}, gridStep, eMinEarth_deg, eMinRelay_deg, ...
    sprintf('Worst at \\phi_{offset}=%d^\\circ  |  Any(OK)=%.2f%%  |  None=%.2f%%', ...
            phiOffsets_deg(idxWorst), pctAny(idxWorst), pctNone(idxWorst)));

function [u_rel, dvec] = buildRing(Nsat, phiOffset_deg, alpha_deg, dL2_km)
    d2r = pi/180;
    u0 = [-1;0;0]; b1 = [0;0;1]; b2 = cross(u0,b1); b2 = b2/norm(b2);
    phi_deg = phiOffset_deg + (0:Nsat-1)*(360/Nsat);
    phi = phi_deg*d2r; alpha = alpha_deg*d2r;
    dvec = dL2_km(:).';

    u_rel = zeros(3,Nsat);
    for k = 1:Nsat
        uk = cos(alpha)*u0 + sin(alpha)*(cos(phi(k))*b1 + sin(phi(k))*b2);
        u_rel(:,k) = uk / norm(uk);
    end
end

function [res, eRelayMax, minRhoVis] = coverage_and_worst(Rpos, rx, ry, rz, nx, ny, nz, eMinRelay_deg, EarthOK, latG)
    [nr, nc] = size(latG);
    Nsat = size(Rpos,2);
    RelayOK_any = false(nr,nc);
    eRelayMax   = -inf(nr,nc);
    minRhoVis   = inf(nr,nc);

    for k = 1:Nsat
        rxR = Rpos(1,k); ryR = Rpos(2,k); rzR = Rpos(3,k); %zenith
        rhox = rxR - rx;  rhoy = ryR - ry;  rhoz = rzR - rz; %los
        rho  = sqrt(rhox.^2 + rhoy.^2 + rhoz.^2); %norm

        ux = rhox ./ rho;  uy = rhoy ./ rho;  uz = rhoz ./ rho; %unit vector of los
        eR = asind( max(min(ux.*nx + uy.*ny + uz.*nz,1),-1) ); % elevation

        vis = (eR >= eMinRelay_deg);
        RelayOK_any = RelayOK_any | vis; % mark visible
        eRelayMax   = max(eRelayMax, eR); % lets update the maximum elevation for each

        rtmp = rho;
        rtmp(~vis) = inf; % if not visble then rule out
        minRhoVis = min(minRhoVis, rtmp); %min over visibles
    end
    minRhoVis(~RelayOK_any) = NaN; % we dont need unvisibles

    EarthOnly =  EarthOK & ~RelayOK_any;
    RelayOnly = ~EarthOK &  RelayOK_any;
    BothOK    =  EarthOK &  RelayOK_any;
    None      = ~EarthOK & ~RelayOK_any;
    AnyOK     = EarthOnly | RelayOnly | BothOK;

    w = cosd(latG); % lets weight each cell
    w = w / sum(w(:));
    res.pct_None  = 100*sum(w(None));
    res.pct_Earth = 100*sum(w(EarthOnly));
    res.pct_Relay = 100*sum(w(RelayOnly));
    res.pct_Both  = 100*sum(w(BothOK));
    res.pct_Any   = 100*sum(w(AnyOK));

    catMap = zeros(nr,nc,'uint8');
    catMap(EarthOnly)=1; catMap(RelayOnly)=2; catMap(BothOK)=3;
    res.catMap = catMap;
end

function m = maxfinite(v)
    v = v(isfinite(v)); % ignore points that have not visible property
    if isempty(v), m = NaN; else, m = max(v); end
end

function m = minfinite(v)
    v = v(isfinite(v));
    if isempty(v), m = NaN; else, m = min(v); end
end

function plotCategoryMap(catMap, gridStep, ~, ~, ~)
    % Grid
    lon = 0:gridStep:360-gridStep;
    lat = -90:gridStep:90;

    % Colors
    cmap = [0.92 0.92 0.92;   % 0 None
            0.20 0.60 1.00;   % 1 Earth only
            1.00 0.55 0.10;   % 2 Relay only
            0.20 0.80 0.25];  % 3 Both
    labels = {'No Coverage','Earth only','Relay only','Both'};

    f = figure('Color','w','Name','Coverage map');
    ax = axes('Parent',f, 'Position',[0.08 0.10 0.73 0.82]);  % [L B W H]
    imagesc(ax, lon, lat, double(catMap));
    set(ax,'YDir','normal'); axis(ax,[0 360 -90 90]); axis(ax,'tight'); box(ax,'on');
    colormap(ax, cmap); caxis(ax,[-0.5 3.5]);
    xlabel(ax,'Longitude (deg)'); ylabel(ax,'Latitude (deg)');
    title(ax,'Coverage map');

    tileW = 0.035; tileH = 0.045;    % figure units
    x0 = 0.85; yTop = 0.86; vgap = 0.035;
    for i = 1:4
        y = yTop - (i-1)*(tileH+vgap);
        annotation(f,'rectangle',[x0 y tileW tileH], ...
            'FaceColor',cmap(i,:), 'EdgeColor','k', 'LineWidth',0.75);
        annotation(f,'textbox',[x0+tileW+0.01 y tileW*3 tileH], ...
            'String',labels{i}, 'EdgeColor','none', ...
            'VerticalAlignment','middle','FontSize',10);
    end
end