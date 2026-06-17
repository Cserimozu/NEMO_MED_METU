% Multi-year event-based offline LIW particle tracking without vertical velocity.
%
% Current parametrization:
%   S >= 39.12
%   14.5 <= T <= 18.0 degC
%   28.60 <= sigma0 <= 29.00 kg m-3
%   120 <= z <= 400 m
%   initialization at the 5-day record closest to 1 March
%   daily offline substepping within 5-day mean velocity fields
%   4-year forward tracking window
%
% Important interpretation:
%   Particles are initialized from LIW-class water, then advected passively
%   at fixed release depth. LIW criteria are used only at initialization.
%   Particles are not removed if their later T/S/sigma0 environment leaves
%   the LIW range.
%
% With model output ending in 2022, a 4-year tracking window allows release
% years through 2018 while keeping an equal forward-tracking duration.

clear; clc;

%% =========================
% CONFIGURATION
% =========================

cfg.rootDir = 'E:\Results10';

% Full-period annual 4-year experiment.
% With 4-year tracking and data ending in 2022, 2018 is the latest valid release year.
cfg.releaseYears = 1992:2018;

cfg.trackYears = 4;
cfg.requireFullTrackingWindow = true;

cfg.runName = 'event_tracking_1992_2018_annual_4yr_dt1_Smin3912_Tmax18_Mar01_no29E_blockedDry60d';

% Output timing.
% Velocity fields are 5-day means, but particles are moved in 1-day substeps.
cfg.dtOutputDays = 5;
cfg.dtSubDays    = 1;

% Release timing.
% Use the actual NetCDF time vector and select the available 5-day record
% closest to 1 March. This keeps seeding in the late-winter / early-spring
% LIW detection phase, rather than relying on a hard-coded record number.
% Change [3 1] to [2 28] for late February or [3 15] for mid-March.
cfg.releaseTargetMonthDay = [3 1];

% Fallback only if the NetCDF time variable cannot be read.
% For 5-day output, record 13 is usually early March.
cfg.releaseFallbackRecordInYear = 13;

% Particle number per source region per release year.
% 27 release years x 3 regions x 200 particles = 16200 particles.
cfg.maxParticlesPerRegion = 200;

% Random seed.
cfg.randomSeed = 42;

% LIW seeding criteria.
cfg.seed.zMin = 120;
cfg.seed.zMax = 400;

% Stricter salinity threshold.
cfg.seed.SMin = 39.12;

cfg.seed.TMin = 14.5;
cfg.seed.TMax = 18.0;

% Density filter.
% Uses TEOS-10 / GSW if available.
cfg.seed.useSigma0 = true;
cfg.seed.sigMin = 28.60;
cfg.seed.sigMax = 29.00;

% Source-region boxes.
cfg.regions = struct([]);

cfg.regions(1).name   = 'Rhodes_Gyre';
cfg.regions(1).lonMin = 27.8;
cfg.regions(1).lonMax = 30.5;
cfg.regions(1).latMin = 34.5;
cfg.regions(1).latMax = 36.8;

cfg.regions(2).name   = 'Antalya_Basin';
cfg.regions(2).lonMin = 29.5;
cfg.regions(2).lonMax = 32.5;
cfg.regions(2).latMin = 35.5;
cfg.regions(2).latMax = 37.0;

cfg.regions(3).name   = 'Cilician_Basin';
cfg.regions(3).lonMin = 32.0;
cfg.regions(3).lonMax = 35.5;
cfg.regions(3).latMin = 35.0;
cfg.regions(3).latMax = 36.8;

% Land/sea mask plotting.
% Derived convention:
%   0 = land/dry
%   1 = sea/wet
cfg.mask.plotLand = true;
cfg.mask.landValue = 0;

% Western open boundary.
% The physical open boundary is the western model edge.
cfg.openBoundary.rimWidth = 10;

% Dry/topographic blocking threshold.
% With cfg.dtSubDays = 1, this means 60 consecutive blocked 1-day steps.
cfg.dryBlock.maxConsecutiveBlockedDays = 60;

% Output folder.
cfg.outDir = fullfile(cfg.rootDir, 'LIW_residence_noW_EVENT_output');

if ~exist(cfg.outDir, 'dir')
    mkdir(cfg.outDir);
end

rng(cfg.randomSeed);

%% =========================
% BUILD YEARLY RECORD LIST
% =========================

fprintf('Building yearly file list...\n');

extraYearsNeeded = ceil(cfg.trackYears);
allYearsNeeded = min(cfg.releaseYears) : max(cfg.releaseYears) + extraYearsNeeded;

records = struct([]);

for yy = allYearsNeeded

    folder = fullfile(cfg.rootDir, num2str(yy));

    Tfile = fullfile(folder, sprintf('NEMO_T_%04d_NaN.nc', yy));
    Ufile = fullfile(folder, sprintf('NEMO_U_%04d_NaN.nc', yy));
    Vfile = fullfile(folder, sprintf('NEMO_V_%04d_NaN.nc', yy));

    if ~isfile(Tfile) || ~isfile(Ufile) || ~isfile(Vfile)
        fprintf('Skipping %d because one or more NEMO files are missing.\n', yy);
        continue;
    end

    if ~exist('varNames', 'var')
        varNames = detectVariableNames(Tfile, Ufile, Vfile);
    end

    nt = getTimeLength(Tfile, varNames.T);

    timeVec = readTimeVector(Tfile, varNames.timeT, nt);

    for tt = 1:nt
        k = numel(records) + 1;
        records(k).year  = yy;
        records(k).tidx  = tt;
        records(k).date  = timeVec(tt);
        records(k).Tfile = Tfile;
        records(k).Ufile = Ufile;
        records(k).Vfile = Vfile;
    end
end

if isempty(records)
    error('No valid yearly files found. Check cfg.rootDir and file names.');
end

fprintf('Total available records: %d\n', numel(records));

fprintf('\nConfigured annual initialization schedule:\n');

for yy = cfg.releaseYears

    [releaseIndex, releaseLabel] = selectReleaseIndex(records, yy, cfg);

    if isempty(releaseIndex)
        fprintf('%d : no release record found\n', yy);
    else
        if isfield(records, 'date') && ~isnat(records(releaseIndex).date)
            fprintf('%d : %s | record %d | %s\n', ...
                yy, releaseLabel, records(releaseIndex).tidx, ...
                datestr(records(releaseIndex).date));
        else
            fprintf('%d : %s | record %d | date unavailable from NetCDF time variable\n', ...
                yy, releaseLabel, records(releaseIndex).tidx);
        end
    end
end

%% =========================
% READ GRIDS
% =========================

fprintf('\nReading grids...\n');

firstT = records(1).Tfile;
firstU = records(1).Ufile;
firstV = records(1).Vfile;

gridT = readGrid(firstT, varNames.lonT, varNames.latT, varNames.depthT);
gridU = readGrid(firstU, varNames.lonU, varNames.latU, varNames.depthU);
gridV = readGrid(firstV, varNames.lonV, varNames.latV, varNames.depthV);

fprintf('T grid: nx=%d, ny=%d, nz=%d\n', ...
    numel(gridT.lon), numel(gridT.lat), numel(gridT.z));

fprintf('U grid: nx=%d, ny=%d, nz=%d\n', ...
    numel(gridU.lon), numel(gridU.lat), numel(gridU.z));

fprintf('V grid: nx=%d, ny=%d, nz=%d\n', ...
    numel(gridV.lon), numel(gridV.lat), numel(gridV.z));

%% =========================
% DERIVE LAND/SEA MASK AND WESTERN OPEN BOUNDARY
% =========================

fprintf('Deriving land/sea mask from NEMO_T zero-land structure...\n');

Tmask0 = read4Drecord(records(1).Tfile, varNames.T, records(1).tidx);

[plotMask, wetMask2D, wetMask3D] = deriveMaskFromTemperature(Tmask0);

cfg.landCheck.lon = gridT.lon;
cfg.landCheck.lat = gridT.lat;
cfg.landCheck.z   = gridT.z;
cfg.landCheck.wetMask3D = wetMask3D;

fprintf('Land/sea and 3D wet/dry masks derived from temperature field.\n');

fprintf('Deriving western open-boundary line from model western edge...\n');

westBoundary = deriveWesternOpenBoundaryFromGrid( ...
    gridT, gridU, wetMask2D, cfg.openBoundary.rimWidth);

cfg.westBoundary = westBoundary;

fprintf('Western open boundary longitude: %.4f\n', westBoundary.lon0);
fprintf('Western open boundary latitude range: %.4f to %.4f\n', ...
    westBoundary.latMin, westBoundary.latMax);

boundaryFile = fullfile(cfg.outDir, ...
    sprintf('western_open_boundary_coordinates_%s.csv', cfg.runName));

boundaryTable = table(westBoundary.lat(:), westBoundary.lon(:), ...
    'VariableNames', {'lat','western_boundary_lon'});

writetable(boundaryTable, boundaryFile);

fprintf('Saved western boundary coordinates:\n%s\n', boundaryFile);

%% =========================
% MAIN ANNUAL RELEASE LOOP
% =========================

allParticles = table();

for yy = cfg.releaseYears

    [releaseIndex, releaseLabel] = selectReleaseIndex(records, yy, cfg);

    if isempty(releaseIndex)
        fprintf('No release record found for %d. Skipping.\n', yy);
        continue;
    end

    maxStepsNeeded = round(cfg.trackYears * 365.25 / cfg.dtOutputDays);
    requiredFinalIndex = releaseIndex + maxStepsNeeded - 1;

    if cfg.requireFullTrackingWindow && requiredFinalIndex > numel(records)
        fprintf(['Skipping release year %d because the full %.2f-year ', ...
                 'tracking window is not available.\n'], yy, cfg.trackYears);
        continue;
    end

    fprintf('\n======================================\n');
    fprintf('Release year: %d | %s | record in year: %d\n', ...
        yy, releaseLabel, records(releaseIndex).tidx);

    if isfield(records, 'date') && ~isnat(records(releaseIndex).date)
        fprintf('Release date from file: %s\n', datestr(records(releaseIndex).date));
    end

    fprintf('Tracking duration: %.2f years\n', cfg.trackYears);
    fprintf('Offline substep: %.1f days\n', cfg.dtSubDays);
    fprintf('SMin: %.2f\n', cfg.seed.SMin);
    fprintf('T range: %.2f to %.2f degC\n', cfg.seed.TMin, cfg.seed.TMax);
    fprintf('sigma0 range: %.2f to %.2f\n', cfg.seed.sigMin, cfg.seed.sigMax);
    fprintf('======================================\n');

    fprintf('Reading T/S for annual seeding...\n');

    T0 = read4Drecord(records(releaseIndex).Tfile, ...
        varNames.T, records(releaseIndex).tidx);

    S0 = read4Drecord(records(releaseIndex).Tfile, ...
        varNames.S, records(releaseIndex).tidx);

    seeds = seedLIWParticles(T0, S0, gridT, cfg, yy, releaseIndex);

    fprintf('Total particles seeded for %d: %d\n', yy, height(seeds));

    if isempty(seeds)
        continue;
    end

    result = trackParticlesEventsNoW(seeds, releaseIndex, records, ...
        gridU, gridV, cfg, varNames);

    allParticles = [allParticles; result]; %#ok<AGROW>

    partialFile = fullfile(cfg.outDir, ...
        sprintf('LIW_event_particles_until_%04d_%s.csv', yy, cfg.runName));

    writetable(allParticles, partialFile);
end

%% =========================
% SAVE OUTPUTS
% =========================

if isempty(allParticles)
    error('No particles were produced.');
end

particleFile = fullfile(cfg.outDir, ...
    sprintf('LIW_event_particles_all_%s.csv', cfg.runName));

yearRegionSummaryFile = fullfile(cfg.outDir, ...
    sprintf('LIW_event_summary_BY_RELEASE_YEAR_AND_REGION_%s.csv', cfg.runName));

yearBasinSummaryFile = fullfile(cfg.outDir, ...
    sprintf('LIW_event_summary_BY_RELEASE_YEAR_BASIN_WIDE_%s.csv', cfg.runName));

regionalSummaryFile = fullfile(cfg.outDir, ...
    sprintf('LIW_event_summary_ALL_YEARS_BY_REGION_%s.csv', cfg.runName));

basinSummaryFile = fullfile(cfg.outDir, ...
    sprintf('LIW_event_summary_ALL_YEARS_BASIN_WIDE_%s.csv', cfg.runName));

writetable(allParticles, particleFile);

[yearRegionSummaryTable, yearBasinSummaryTable, ...
 regionalSummaryTable, basinSummaryTable] = summarizeEventsMultiYear(allParticles);

writetable(yearRegionSummaryTable, yearRegionSummaryFile);
writetable(yearBasinSummaryTable, yearBasinSummaryFile);
writetable(regionalSummaryTable, regionalSummaryFile);
writetable(basinSummaryTable, basinSummaryFile);

fprintf('\nSaved particle table:\n%s\n', particleFile);
fprintf('Saved yearly regional summary table:\n%s\n', yearRegionSummaryFile);
fprintf('Saved yearly basin-wide summary table:\n%s\n', yearBasinSummaryFile);
fprintf('Saved all-years regional summary table:\n%s\n', regionalSummaryFile);
fprintf('Saved all-years basin-wide summary table:\n%s\n', basinSummaryFile);

fprintf('\nYearly regional summary:\n');
disp(yearRegionSummaryTable);

fprintf('\nAll-years regional summary:\n');
disp(regionalSummaryTable);

fprintf('\nAll-years basin-wide summary:\n');
disp(basinSummaryTable);

%% =========================
% FIGURES
% =========================

plotEventMap(allParticles, gridT, plotMask, cfg, ...
    'source_exit_months', ...
    'Source-region exit time [months]', ...
    fullfile(cfg.outDir, sprintf('map_source_exit_%s.png', cfg.runName)));

plotEventMap(allParticles, gridT, plotMask, cfg, ...
    'west_boundary_exit_months', ...
    'Western open-boundary exit time [months]', ...
    fullfile(cfg.outDir, sprintf('map_west_boundary_exit_%s.png', cfg.runName)));

plotEventHistogram(allParticles.source_exit_months, ...
    'Source-region exit time [months]', ...
    fullfile(cfg.outDir, sprintf('hist_source_exit_%s.png', cfg.runName)));

plotEventHistogram(allParticles.west_boundary_exit_months, ...
    'Western open-boundary exit time [months]', ...
    fullfile(cfg.outDir, sprintf('hist_west_boundary_exit_%s.png', cfg.runName)));

plotEventHistogram(allParticles.max_dry_blocked_days, ...
    'Maximum consecutive dry-blocked time [days]', ...
    fullfile(cfg.outDir, sprintf('hist_max_dry_blocked_days_%s.png', cfg.runName)));

plotRegionalBoxplot(allParticles, ...
    'source_exit_months', ...
    'Source-region exit time [months]', ...
    fullfile(cfg.outDir, sprintf('box_source_exit_by_region_%s.png', cfg.runName)));

plotRegionalBoxplot(allParticles, ...
    'west_boundary_exit_months', ...
    'Western open-boundary exit time [months]', ...
    fullfile(cfg.outDir, sprintf('box_west_boundary_exit_by_region_%s.png', cfg.runName)));

plotRegionalBoxplot(allParticles, ...
    'max_dry_blocked_days', ...
    'Maximum consecutive dry-blocked time [days]', ...
    fullfile(cfg.outDir, sprintf('box_max_dry_blocked_by_region_%s.png', cfg.runName)));

%% =========================
% RELEASE-YEAR TIME SERIES
% =========================

plotAnnualMetricByRegion(yearRegionSummaryTable, ...
    'median_source_exit_days', ...
    'Median source-region exit time [days]', ...
    fullfile(cfg.outDir, sprintf('annual_median_source_exit_by_region_%s.png', cfg.runName)));

plotAnnualMetricByRegion(yearRegionSummaryTable, ...
    'median_west_boundary_exit_days', ...
    'Median western-boundary exit time [days]', ...
    fullfile(cfg.outDir, sprintf('annual_median_west_exit_by_region_%s.png', cfg.runName)));

plotAnnualMetricByRegion(yearRegionSummaryTable, ...
    'pct_west_boundary_exit', ...
    'Particles exiting western boundary [%]', ...
    fullfile(cfg.outDir, sprintf('annual_pct_west_exit_by_region_%s.png', cfg.runName)));

plotAnnualMetricByRegion(yearRegionSummaryTable, ...
    'pct_terminal_max_time', ...
    'Particles retained/censored at end [%]', ...
    fullfile(cfg.outDir, sprintf('annual_pct_retained_by_region_%s.png', cfg.runName)));

plotAnnualMetricByRegion(yearRegionSummaryTable, ...
    'pct_terminal_topographic_trap', ...
    'Coastal/topographic trap [%]', ...
    fullfile(cfg.outDir, sprintf('annual_pct_topographic_trap_by_region_%s.png', cfg.runName)));

fprintf('\nDone.\n');

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function [releaseIndex, releaseLabel] = selectReleaseIndex(records, yy, cfg)

    releaseIndex = [];
    releaseLabel = '';

    idxYear = find([records.year] == yy);

    if isempty(idxYear)
        return;
    end

    % Prefer the actual NetCDF dates, because the meaning of record numbers
    % can shift slightly depending on time-coordinate convention.
    if isfield(records, 'date')
        dateVec = [records(idxYear).date];
        validDate = ~isnat(dateVec);

        if any(validDate)
            idxValid = idxYear(validDate);
            dateVecValid = [records(idxValid).date];

            targetDate = datetime(yy, cfg.releaseTargetMonthDay(1), ...
                                      cfg.releaseTargetMonthDay(2));

            dtDays = abs(days(dateVecValid - targetDate));
            [~, ii] = min(dtDays);

            releaseIndex = idxValid(ii);
            releaseLabel = sprintf('target %02d-%02d closest record', ...
                cfg.releaseTargetMonthDay(1), cfg.releaseTargetMonthDay(2));
            return;
        end
    end

    % Fallback for files without a readable time variable.
    if isfield(cfg, 'releaseFallbackRecordInYear')
        releaseIndex = find([records.year] == yy & ...
                            [records.tidx] == cfg.releaseFallbackRecordInYear, 1);
        releaseLabel = sprintf('fallback record %d', ...
            cfg.releaseFallbackRecordInYear);
    end
end


function varNames = detectVariableNames(Tfile, Ufile, Vfile)

    varNames.lonT   = pickVar(Tfile, {'nav_lon','lon','longitude'});
    varNames.latT   = pickVar(Tfile, {'nav_lat','lat','latitude'});
    varNames.depthT = pickVar(Tfile, {'deptht','depth','depth_0','gdept_0'});

    varNames.lonU   = pickVar(Ufile, {'nav_lon','lon','longitude'});
    varNames.latU   = pickVar(Ufile, {'nav_lat','lat','latitude'});
    varNames.depthU = pickVar(Ufile, {'depthu','deptht','depth','depth_0','gdept_0'});

    varNames.lonV   = pickVar(Vfile, {'nav_lon','lon','longitude'});
    varNames.latV   = pickVar(Vfile, {'nav_lat','lat','latitude'});
    varNames.depthV = pickVar(Vfile, {'depthv','deptht','depth','depth_0','gdept_0'});

    varNames.T = pickVar(Tfile, {'votemper','thetao','temp','temperature','toce'});
    varNames.S = pickVar(Tfile, {'vosaline','so','salt','salinity','soce'});
    varNames.U = pickVar(Ufile, {'vozocrtx','uo','u','u_velocity'});
    varNames.V = pickVar(Vfile, {'vomecrty','vo','v','v_velocity'});

    varNames.timeT = pickOptionalVar(Tfile, {'time_counter','time','time_centered'});

    fprintf('\nDetected variables:\n');
    disp(varNames);
end

function name = pickVar(ncfile, candidates)

    info = ncinfo(ncfile);
    names = string({info.Variables.Name});

    for i = 1:numel(candidates)
        idx = strcmpi(names, candidates{i});
        if any(idx)
            name = char(names(find(idx, 1)));
            return;
        end
    end

    fprintf('\nVariables in file:\n%s\n', ncfile);
    disp(names');

    error('Could not find any candidate variable: %s', ...
        strjoin(candidates, ', '));
end

function name = pickOptionalVar(ncfile, candidates)

    info = ncinfo(ncfile);
    names = string({info.Variables.Name});

    name = '';

    for i = 1:numel(candidates)
        idx = strcmpi(names, candidates{i});
        if any(idx)
            name = char(names(find(idx, 1)));
            return;
        end
    end
end

function tvec = readTimeVector(ncfile, timeVar, nt)

    tvec = NaT(nt, 1);

    if isempty(timeVar)
        return;
    end

    try
        raw = double(squeeze(ncread(ncfile, timeVar)));
    catch
        return;
    end

    if isempty(raw)
        return;
    end

    raw = raw(:);

    if numel(raw) < nt
        return;
    end

    raw = raw(1:nt);

    try
        units = ncreadatt(ncfile, timeVar, 'units');
    catch
        return;
    end

    if isstring(units)
        units = char(units);
    end

    tok = regexp(units, ...
        '(seconds|second|sec|s|days|day|d|hours|hour|hr|h)\s+since\s+(\d{4}-\d{2}-\d{2})(?:[ T](\d{2}:\d{2}:\d{2}))?', ...
        'tokens', 'once');

    if isempty(tok)
        return;
    end

    unitStr = lower(tok{1});
    dateStr = tok{2};

    if numel(tok) >= 3 && ~isempty(tok{3})
        timeStr = tok{3};
    else
        timeStr = '00:00:00';
    end

    try
        base = datetime([dateStr ' ' timeStr], ...
            'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    catch
        return;
    end

    switch unitStr
        case {'seconds','second','sec','s'}
            tvec = base + seconds(raw);
        case {'hours','hour','hr','h'}
            tvec = base + hours(raw);
        case {'days','day','d'}
            tvec = base + days(raw);
        otherwise
            tvec = NaT(nt, 1);
    end
end

function nt = getTimeLength(ncfile, varName)

    info = ncinfo(ncfile, varName);
    sz = info.Size;

    if numel(sz) == 4
        nt = sz(4);
    else
        nt = 1;
    end
end

function grid = readGrid(ncfile, lonVar, latVar, depthVar)

    lon = double(squeeze(ncread(ncfile, lonVar)));
    lat = double(squeeze(ncread(ncfile, latVar)));
    z   = double(squeeze(ncread(ncfile, depthVar)));

    if ~isvector(lon)
        nx = size(lon,1);
        ny = size(lon,2);

        lonAxis = squeeze(lon(:, round(ny/2)));
        latAxis = squeeze(lat(round(nx/2), :));
    else
        lonAxis = lon(:);
        latAxis = lat(:);
    end

    if ndims(z) == 3
        z = squeeze(z(1,1,:));
    elseif ndims(z) == 4
        z = squeeze(z(1,1,:,1));
    end

    grid.lon = double(lonAxis(:));
    grid.lat = double(latAxis(:));
    grid.z   = double(z(:));

    if any(diff(grid.lon) <= 0)
        error('Longitude axis is not strictly increasing in %s.', ncfile);
    end

    if any(diff(grid.lat) <= 0)
        error('Latitude axis is not strictly increasing in %s.', ncfile);
    end

    if any(diff(grid.z) <= 0)
        error('Depth axis is not strictly increasing in %s.', ncfile);
    end
end

function A = read4Drecord(ncfile, varName, tidx)

    info = ncinfo(ncfile, varName);
    sz = info.Size;

    if numel(sz) == 4
        A = ncread(ncfile, varName, [1 1 1 tidx], [Inf Inf Inf 1]);
        A = squeeze(A);
    elseif numel(sz) == 3
        A = ncread(ncfile, varName);
        A = squeeze(A);
    else
        error('Unsupported dimensions for %s in %s', varName, ncfile);
    end

    A = double(A);

    try
        fv = ncreadatt(ncfile, varName, '_FillValue');
        A(A == double(fv)) = NaN;
    catch
    end

    try
        mv = ncreadatt(ncfile, varName, 'missing_value');
        A(A == double(mv)) = NaN;
    catch
    end

    A(abs(A) > 1e19) = NaN;
end

function [lsm2D, wetMask2D, wetMask3D] = deriveMaskFromTemperature(T)

    tol = 1e-12;

    wetMask3D = isfinite(T) & abs(T) > tol;
    wetMask2D = any(wetMask3D, 3);

    lsm2D = double(wetMask2D);
end

function westBoundary = deriveWesternOpenBoundaryFromGrid(gridT, gridU, wetMask2D, rimWidth)

    nx = numel(gridT.lon);
    ny = numel(gridT.lat);

    if ~isequal(size(wetMask2D), [nx ny])
        error('wetMask2D size does not match gridT dimensions.');
    end

    rimWidth = min(rimWidth, nx);

    wetWestRows = any(wetMask2D(1:rimWidth, :), 1);

    if ~any(wetWestRows)
        warning('No wet points found in western rim. Using full latitude range.');
        wbLat = gridT.lat(:);
    else
        wbLat = gridT.lat(wetWestRows);
        wbLat = wbLat(:);
    end

    westBoundary.lon0 = min(gridU.lon);
    westBoundary.lat  = wbLat;
    westBoundary.lon  = westBoundary.lon0 .* ones(size(wbLat));

    westBoundary.latMin = min(wbLat);
    westBoundary.latMax = max(wbLat);
    westBoundary.lonMin = westBoundary.lon0;
    westBoundary.lonMax = westBoundary.lon0;
end

function seeds = seedLIWParticles(T, S, gridT, cfg, releaseYear, releaseIndex)

    [LON, LAT] = ndgrid(gridT.lon, gridT.lat);

    lon0 = [];
    lat0 = [];
    z0 = [];
    sourceID = [];
    sourceName = strings(0,1);

    useSigmaNow = cfg.seed.useSigma0 && exist('gsw_SA_from_SP', 'file') == 2;

    if cfg.seed.useSigma0 && ~useSigmaNow
        warning('GSW toolbox not found. Seeding will use T/S/depth only, without sigma0.');
    end

    for rr = 1:numel(cfg.regions)

        reg = cfg.regions(rr);

        regionMask = LON >= reg.lonMin & LON <= reg.lonMax & ...
                     LAT >= reg.latMin & LAT <= reg.latMax;

        candLon = [];
        candLat = [];
        candZ = [];

        for kk = 1:numel(gridT.z)

            z = gridT.z(kk);

            if z < cfg.seed.zMin || z > cfg.seed.zMax
                continue;
            end

            Tk = T(:,:,kk);
            Sk = S(:,:,kk);

            mask = regionMask & ...
                   isfinite(Tk) & isfinite(Sk) & ...
                   abs(Tk) > 1e-12 & abs(Sk) > 1e-12 & ...
                   Sk >= cfg.seed.SMin & ...
                   Tk >= cfg.seed.TMin & Tk <= cfg.seed.TMax;

            if useSigmaNow
                sig0 = calcSigma0(Sk, Tk, z, LON, LAT);
                mask = mask & sig0 >= cfg.seed.sigMin & ...
                              sig0 <= cfg.seed.sigMax;
            end

            [ii, jj] = find(mask);

            if ~isempty(ii)
                candLon = [candLon; gridT.lon(ii)]; %#ok<AGROW>
                candLat = [candLat; gridT.lat(jj)]; %#ok<AGROW>
                candZ   = [candZ; repmat(z, numel(ii), 1)]; %#ok<AGROW>
            end
        end

        nCand = numel(candLon);

        if nCand == 0
            fprintf('  %-18s : no LIW seeds\n', reg.name);
            continue;
        end

        if nCand > cfg.maxParticlesPerRegion
            pick = randperm(nCand, cfg.maxParticlesPerRegion);
        else
            pick = 1:nCand;
        end

        candLon = candLon(pick);
        candLat = candLat(pick);
        candZ   = candZ(pick);

        nPick = numel(candLon);

        lon0 = [lon0; candLon]; %#ok<AGROW>
        lat0 = [lat0; candLat]; %#ok<AGROW>
        z0   = [z0; candZ]; %#ok<AGROW>
        sourceID = [sourceID; repmat(rr, nPick, 1)]; %#ok<AGROW>
        sourceName = [sourceName; repmat(string(reg.name), nPick, 1)]; %#ok<AGROW>

        fprintf('  %-18s : %d seeds from %d candidates\n', ...
            reg.name, nPick, nCand);
    end

    n = numel(lon0);

    seeds = table();
    seeds.release_year  = repmat(releaseYear, n, 1);
    seeds.release_index = repmat(releaseIndex, n, 1);
    seeds.source_id     = sourceID;
    seeds.source_region = sourceName;
    seeds.lon0 = lon0;
    seeds.lat0 = lat0;
    seeds.z0   = z0;
end

function result = trackParticlesEventsNoW(seeds, releaseIndex, records, gridU, gridV, cfg, varNames)

    N = height(seeds);

    lon = seeds.lon0;
    lat = seeds.lat0;
    z   = seeds.z0;

    alive = true(N,1);
    ageDays = zeros(N,1);

    sourceExitDays = NaN(N,1);
    sourceExitLon  = NaN(N,1);
    sourceExitLat  = NaN(N,1);

    westBoundaryExitDays = NaN(N,1);
    westBoundaryExitLon  = NaN(N,1);
    westBoundaryExitLat  = NaN(N,1);

    terminalDays = NaN(N,1);
    terminalReason = strings(N,1);
    terminalLon = NaN(N,1);
    terminalLat = NaN(N,1);
    terminalZ   = NaN(N,1);

    dryBlockedDays = zeros(N,1);
    maxDryBlockedDays = zeros(N,1);
    nDryBlockedAttempts = zeros(N,1);

    R = 6371000;
    dtSec = cfg.dtSubDays * 86400;
    nSub = max(1, round(cfg.dtOutputDays / cfg.dtSubDays));

    maxSteps = round(cfg.trackYears * 365.25 / cfg.dtOutputDays);
    finalIndex = min(numel(records), releaseIndex + maxSteps - 1);

    if isfield(cfg, 'landCheck') && ~isempty(cfg.landCheck.wetMask3D)
        Fwet = griddedInterpolant( ...
            {cfg.landCheck.lon, cfg.landCheck.lat, cfg.landCheck.z}, ...
            double(cfg.landCheck.wetMask3D), ...
            'nearest', 'none');
    else
        Fwet = [];
    end

    for recIndex = releaseIndex:finalIndex

        fprintf('  Advecting record %d/%d | year %d step %d | alive %d\n', ...
            recIndex, finalIndex, records(recIndex).year, ...
            records(recIndex).tidx, sum(alive));

        if ~any(alive)
            break;
        end

        U = read4Drecord(records(recIndex).Ufile, ...
            varNames.U, records(recIndex).tidx);

        V = read4Drecord(records(recIndex).Vfile, ...
            varNames.V, records(recIndex).tidx);

        FU = griddedInterpolant({gridU.lon, gridU.lat, gridU.z}, ...
            U, 'linear', 'none');

        FV = griddedInterpolant({gridV.lon, gridV.lat, gridV.z}, ...
            V, 'linear', 'none');

        for sub = 1:nSub

            idx = find(alive);

            if isempty(idx)
                break;
            end

            oldLon = lon;
            oldLat = lat;

            ui = FU(lon(idx), lat(idx), z(idx));
            vi = FV(lon(idx), lat(idx), z(idx));

            bad = ~isfinite(ui) | ~isfinite(vi);

            if any(bad)
                deadIdx = idx(bad);
                [alive, terminalDays, terminalReason, terminalLon, terminalLat, terminalZ] = ...
                    markTerminal(deadIdx, alive, terminalDays, terminalReason, ...
                    terminalLon, terminalLat, terminalZ, ageDays, lon, lat, z, ...
                    "land_or_missing_velocity");
            end

            goodIdx = idx(~bad);

            if isempty(goodIdx)
                continue;
            end

            ui = ui(~bad);
            vi = vi(~bad);

            latRad = deg2rad(lat(goodIdx));

            dlon = rad2deg((ui .* dtSec) ./ (R .* cos(latRad)));
            dlat = rad2deg((vi .* dtSec) ./ R);

            lon(goodIdx) = lon(goodIdx) + dlon;
            lat(goodIdx) = lat(goodIdx) + dlat;

            ageDays(goodIdx) = ageDays(goodIdx) + cfg.dtSubDays;

            westOpenBoundaryExit = detectWesternOpenBoundaryExit( ...
                oldLon, oldLat, lon, lat, alive, cfg.westBoundary);

            if any(westOpenBoundaryExit)

                westBoundaryExitDays(westOpenBoundaryExit) = ageDays(westOpenBoundaryExit);
                westBoundaryExitLon(westOpenBoundaryExit)  = lon(westOpenBoundaryExit);
                westBoundaryExitLat(westOpenBoundaryExit)  = lat(westOpenBoundaryExit);

                for rr = 1:numel(cfg.regions)

                    reg = cfg.regions(rr);

                    idxReg = westOpenBoundaryExit & ...
                             seeds.source_id == rr & ...
                             isnan(sourceExitDays);

                    if ~any(idxReg)
                        continue;
                    end

                    inside = lon >= reg.lonMin & lon <= reg.lonMax & ...
                             lat >= reg.latMin & lat <= reg.latMax;

                    outRegion = idxReg & ~inside;

                    if any(outRegion)
                        sourceExitDays(outRegion) = ageDays(outRegion);
                        sourceExitLon(outRegion)  = lon(outRegion);
                        sourceExitLat(outRegion)  = lat(outRegion);
                    end
                end

                deadIdx = find(westOpenBoundaryExit);

                [alive, terminalDays, terminalReason, terminalLon, terminalLat, terminalZ] = ...
                    markTerminal(deadIdx, alive, terminalDays, terminalReason, ...
                    terminalLon, terminalLat, terminalZ, ageDays, lon, lat, z, ...
                    "west_open_boundary_exit");
            end

            nonWesternDomainError = alive & ...
                (lon > max(gridU.lon) | ...
                 lat < min(gridU.lat) | ...
                 lat > max(gridU.lat));

            if any(nonWesternDomainError)
                deadIdx = find(nonWesternDomainError);

                [alive, terminalDays, terminalReason, terminalLon, terminalLat, terminalZ] = ...
                    markTerminal(deadIdx, alive, terminalDays, terminalReason, ...
                    terminalLon, terminalLat, terminalZ, ageDays, lon, lat, z, ...
                    "nonwestern_domain_error");
            end

            if ~isempty(Fwet)

                idxAlive = find(alive);

                if ~isempty(idxAlive)
                    wetHere = Fwet(lon(idxAlive), lat(idxAlive), z(idxAlive));
                    hitDry = ~isfinite(wetHere) | wetHere < 0.5;

                    if any(hitDry)

                        blockedIdx = idxAlive(hitDry);

                        lon(blockedIdx) = oldLon(blockedIdx);
                        lat(blockedIdx) = oldLat(blockedIdx);

                        dryBlockedDays(blockedIdx) = ...
                            dryBlockedDays(blockedIdx) + cfg.dtSubDays;

                        maxDryBlockedDays(blockedIdx) = max( ...
                            maxDryBlockedDays(blockedIdx), ...
                            dryBlockedDays(blockedIdx));

                        nDryBlockedAttempts(blockedIdx) = ...
                            nDryBlockedAttempts(blockedIdx) + 1;
                    end

                    if any(~hitDry)
                        movedWetIdx = idxAlive(~hitDry);
                        dryBlockedDays(movedWetIdx) = 0;
                    end
                end
            end

            topographicTrap = alive & ...
                dryBlockedDays >= cfg.dryBlock.maxConsecutiveBlockedDays;

            if any(topographicTrap)
                deadIdx = find(topographicTrap);

                [alive, terminalDays, terminalReason, terminalLon, terminalLat, terminalZ] = ...
                    markTerminal(deadIdx, alive, terminalDays, terminalReason, ...
                    terminalLon, terminalLat, terminalZ, ageDays, lon, lat, z, ...
                    "coastal_topographic_trap");
            end

            for rr = 1:numel(cfg.regions)

                reg = cfg.regions(rr);

                idxReg = alive & seeds.source_id == rr & isnan(sourceExitDays);

                if ~any(idxReg)
                    continue;
                end

                inside = lon >= reg.lonMin & lon <= reg.lonMax & ...
                         lat >= reg.latMin & lat <= reg.latMax;

                outRegion = idxReg & ~inside;

                if any(outRegion)
                    sourceExitDays(outRegion) = ageDays(outRegion);
                    sourceExitLon(outRegion)  = lon(outRegion);
                    sourceExitLat(outRegion)  = lat(outRegion);
                end
            end
        end
    end

    stillAlive = find(alive);

    if ~isempty(stillAlive)
        [alive, terminalDays, terminalReason, terminalLon, terminalLat, terminalZ] = ...
            markTerminal(stillAlive, alive, terminalDays, terminalReason, ...
            terminalLon, terminalLat, terminalZ, ageDays, lon, lat, z, ...
            "max_tracking_time");
    end

    result = seeds;

    result.source_exit_days   = sourceExitDays;
    result.source_exit_months = sourceExitDays ./ 30.4375;
    result.source_exit_lon    = sourceExitLon;
    result.source_exit_lat    = sourceExitLat;

    result.west_boundary_exit_days   = westBoundaryExitDays;
    result.west_boundary_exit_months = westBoundaryExitDays ./ 30.4375;
    result.west_boundary_exit_lon    = westBoundaryExitLon;
    result.west_boundary_exit_lat    = westBoundaryExitLat;

    result.terminal_days   = terminalDays;
    result.terminal_months = terminalDays ./ 30.4375;
    result.terminal_reason = terminalReason;
    result.terminal_lon    = terminalLon;
    result.terminal_lat    = terminalLat;
    result.terminal_z      = terminalZ;

    result.n_dry_blocked_attempts = nDryBlockedAttempts;
    result.max_dry_blocked_days   = maxDryBlockedDays;
end

function crossed = detectWesternOpenBoundaryExit(oldLon, oldLat, newLon, newLat, activeMask, westBoundary)

    latMid = 0.5 .* (oldLat + newLat);

    inBoundaryLatRange = latMid >= westBoundary.latMin & ...
                         latMid <= westBoundary.latMax;

    crossed = activeMask & ...
              inBoundaryLatRange & ...
              oldLon >= westBoundary.lon0 & ...
              newLon <  westBoundary.lon0;
end

function [alive, terminalDays, terminalReason, terminalLon, terminalLat, terminalZ] = ...
    markTerminal(deadIdx, alive, terminalDays, terminalReason, ...
    terminalLon, terminalLat, terminalZ, ageDays, lon, lat, z, reason)

    deadIdx = deadIdx(alive(deadIdx));

    if isempty(deadIdx)
        return;
    end

    alive(deadIdx) = false;
    terminalDays(deadIdx) = ageDays(deadIdx);
    terminalReason(deadIdx) = reason;
    terminalLon(deadIdx) = lon(deadIdx);
    terminalLat(deadIdx) = lat(deadIdx);
    terminalZ(deadIdx)   = z(deadIdx);
end

function sig0 = calcSigma0(SP, T, z, LON, LAT)

    p  = gsw_p_from_z(-z .* ones(size(SP)), LAT);
    SA = gsw_SA_from_SP(SP, p, LON, LAT);
    CT = gsw_CT_from_t(SA, T, p);
    sig0 = gsw_sigma0(SA, CT);
end

function [yearRegionSummaryTable, yearBasinSummaryTable, ...
          regionalSummaryTable, basinSummaryTable] = summarizeEventsMultiYear(P)

    years = unique(P.release_year);
    regions = unique(P.source_region);

    yearRegionSummaryTable = table();
    yearBasinSummaryTable = table();
    regionalSummaryTable = table();

    for yy = years(:)'

        for rr = 1:numel(regions)

            r = regions(rr);
            idx = P.release_year == yy & P.source_region == r;

            if ~any(idx)
                continue;
            end

            T = summarizeOneGroup(P, idx, r);
            T.release_year = yy;
            T = movevars(T, 'release_year', 'Before', 'group');

            yearRegionSummaryTable = [yearRegionSummaryTable; T]; %#ok<AGROW>
        end

        idxYear = P.release_year == yy;

        if any(idxYear)
            T = summarizeOneGroup(P, idxYear, "BASIN_WIDE");
            T.release_year = yy;
            T = movevars(T, 'release_year', 'Before', 'group');

            yearBasinSummaryTable = [yearBasinSummaryTable; T]; %#ok<AGROW>
        end
    end

    for rr = 1:numel(regions)

        r = regions(rr);
        idx = P.source_region == r;

        T = summarizeOneGroup(P, idx, r);
        T.release_year = NaN;
        T = movevars(T, 'release_year', 'Before', 'group');

        regionalSummaryTable = [regionalSummaryTable; T]; %#ok<AGROW>
    end

    basinSummaryTable = summarizeOneGroup(P, true(height(P),1), "BASIN_WIDE");
    basinSummaryTable.release_year = NaN;
    basinSummaryTable = movevars(basinSummaryTable, 'release_year', 'Before', 'group');
end

function T = summarizeOneGroup(P, idx, groupName)

    sx = P.source_exit_days(idx);
    wx = P.west_boundary_exit_days(idx);
    term = P.terminal_reason(idx);

    n_particles = sum(idx);

    n_source_exit = sum(isfinite(sx));
    n_west_boundary_exit = sum(isfinite(wx));

    n_terminal_west_boundary = sum(term == "west_open_boundary_exit");
    n_terminal_topographic_trap = sum(term == "coastal_topographic_trap");
    n_terminal_land_missing = sum(term == "land_or_missing_velocity");
    n_terminal_nonwestern_error = sum(term == "nonwestern_domain_error");
    n_terminal_max_time = sum(term == "max_tracking_time");

    pct_source_exit = 100 * n_source_exit / n_particles;
    pct_west_boundary_exit = 100 * n_west_boundary_exit / n_particles;

    pct_terminal_west_boundary = 100 * n_terminal_west_boundary / n_particles;
    pct_terminal_topographic_trap = 100 * n_terminal_topographic_trap / n_particles;
    pct_terminal_land_missing = 100 * n_terminal_land_missing / n_particles;
    pct_terminal_nonwestern_error = 100 * n_terminal_nonwestern_error / n_particles;
    pct_terminal_max_time = 100 * n_terminal_max_time / n_particles;

    n_no_source_exit_recorded = sum(~isfinite(sx));
    n_no_west_boundary_exit_recorded = sum(~isfinite(wx));

    pct_no_source_exit_recorded = 100 * n_no_source_exit_recorded / n_particles;
    pct_no_west_boundary_exit_recorded = ...
        100 * n_no_west_boundary_exit_recorded / n_particles;

    median_source_exit_days = safeMedian(sx);
    p25_source_exit_days = safePrctile(sx, 25);
    p75_source_exit_days = safePrctile(sx, 75);

    median_source_exit_months = median_source_exit_days / 30.4375;
    p25_source_exit_months = p25_source_exit_days / 30.4375;
    p75_source_exit_months = p75_source_exit_days / 30.4375;

    median_west_boundary_exit_days = safeMedian(wx);
    p25_west_boundary_exit_days = safePrctile(wx, 25);
    p75_west_boundary_exit_days = safePrctile(wx, 75);

    median_west_boundary_exit_months = median_west_boundary_exit_days / 30.4375;
    p25_west_boundary_exit_months = p25_west_boundary_exit_days / 30.4375;
    p75_west_boundary_exit_months = p75_west_boundary_exit_days / 30.4375;

    source_exit_within_3mo_pct = 100 * sum(isfinite(sx) & sx <= 91.3125) / n_particles;
    source_exit_within_6mo_pct = 100 * sum(isfinite(sx) & sx <= 182.625) / n_particles;
    source_exit_within_1yr_pct = 100 * sum(isfinite(sx) & sx <= 365.25) / n_particles;
    source_exit_within_2yr_pct = 100 * sum(isfinite(sx) & sx <= 730.5) / n_particles;
    source_exit_within_4yr_pct = 100 * sum(isfinite(sx) & sx <= 1461.0) / n_particles;

    west_exit_within_3mo_pct = 100 * sum(isfinite(wx) & wx <= 91.3125) / n_particles;
    west_exit_within_6mo_pct = 100 * sum(isfinite(wx) & wx <= 182.625) / n_particles;
    west_exit_within_1yr_pct = 100 * sum(isfinite(wx) & wx <= 365.25) / n_particles;
    west_exit_within_2yr_pct = 100 * sum(isfinite(wx) & wx <= 730.5) / n_particles;
    west_exit_within_4yr_pct = 100 * sum(isfinite(wx) & wx <= 1461.0) / n_particles;
    blockedAttempts = P.n_dry_blocked_attempts(idx);
    maxBlockedDays  = P.max_dry_blocked_days(idx);

    pct_ever_dry_blocked = 100 * sum(blockedAttempts > 0) / n_particles;
    median_max_dry_blocked_days = safeMedian(maxBlockedDays);
    p75_max_dry_blocked_days = safePrctile(maxBlockedDays, 75);
    max_max_dry_blocked_days = max(maxBlockedDays);

    T = table();

    T.group = string(groupName);
    T.n_particles = n_particles;

    T.n_source_exit = n_source_exit;
    T.pct_source_exit = pct_source_exit;
    T.median_source_exit_days = median_source_exit_days;
    T.p25_source_exit_days = p25_source_exit_days;
    T.p75_source_exit_days = p75_source_exit_days;
    T.median_source_exit_months = median_source_exit_months;
    T.p25_source_exit_months = p25_source_exit_months;
    T.p75_source_exit_months = p75_source_exit_months;
    T.source_exit_within_3mo_pct = source_exit_within_3mo_pct;
    T.source_exit_within_6mo_pct = source_exit_within_6mo_pct;
    T.source_exit_within_1yr_pct = source_exit_within_1yr_pct;
    T.source_exit_within_2yr_pct = source_exit_within_2yr_pct;
    T.source_exit_within_4yr_pct = source_exit_within_4yr_pct;

    T.n_west_boundary_exit = n_west_boundary_exit;
    T.pct_west_boundary_exit = pct_west_boundary_exit;
    T.median_west_boundary_exit_days = median_west_boundary_exit_days;
    T.p25_west_boundary_exit_days = p25_west_boundary_exit_days;
    T.p75_west_boundary_exit_days = p75_west_boundary_exit_days;
    T.median_west_boundary_exit_months = median_west_boundary_exit_months;
    T.p25_west_boundary_exit_months = p25_west_boundary_exit_months;
    T.p75_west_boundary_exit_months = p75_west_boundary_exit_months;
    T.west_exit_within_3mo_pct = west_exit_within_3mo_pct;
    T.west_exit_within_6mo_pct = west_exit_within_6mo_pct;
    T.west_exit_within_1yr_pct = west_exit_within_1yr_pct;
    T.west_exit_within_2yr_pct = west_exit_within_2yr_pct;
    T.west_exit_within_4yr_pct = west_exit_within_4yr_pct;
    T.n_terminal_west_boundary = n_terminal_west_boundary;
    T.pct_terminal_west_boundary = pct_terminal_west_boundary;

    T.n_terminal_topographic_trap = n_terminal_topographic_trap;
    T.pct_terminal_topographic_trap = pct_terminal_topographic_trap;

    T.n_terminal_land_missing = n_terminal_land_missing;
    T.pct_terminal_land_missing = pct_terminal_land_missing;

    T.n_terminal_nonwestern_error = n_terminal_nonwestern_error;
    T.pct_terminal_nonwestern_error = pct_terminal_nonwestern_error;

    T.n_terminal_max_time = n_terminal_max_time;
    T.pct_terminal_max_time = pct_terminal_max_time;

    T.n_no_source_exit_recorded = n_no_source_exit_recorded;
    T.pct_no_source_exit_recorded = pct_no_source_exit_recorded;

    T.n_no_west_boundary_exit_recorded = n_no_west_boundary_exit_recorded;
    T.pct_no_west_boundary_exit_recorded = pct_no_west_boundary_exit_recorded;

    T.pct_ever_dry_blocked = pct_ever_dry_blocked;
    T.median_max_dry_blocked_days = median_max_dry_blocked_days;
    T.p75_max_dry_blocked_days = p75_max_dry_blocked_days;
    T.max_max_dry_blocked_days = max_max_dry_blocked_days;
end

function m = safeMedian(x)
    x = x(isfinite(x));
    if isempty(x)
        m = NaN;
    else
        m = median(x);
    end
end

function q = safePrctile(x, p)
    x = x(isfinite(x));
    if isempty(x)
        q = NaN;
    else
        q = prctile(x, p);
    end
end

function plotEventMap(P, gridT, lsm, cfg, valueName, cbLabel, outFile)

    fig = figure('Color', 'w');
    hold on;

    [LON, LAT] = ndgrid(gridT.lon, gridT.lat);

    if cfg.mask.plotLand && ~isempty(lsm)
        landMask = lsm == cfg.mask.landValue;
        scatter(LON(landMask), LAT(landMask), 4, [0.75 0.75 0.75], 'filled');
    end

    if isfield(cfg, 'westBoundary')
        plot(cfg.westBoundary.lon, cfg.westBoundary.lat, 'k-', 'LineWidth', 1.5);
    end

    vals = P.(valueName);

    hasEvent = isfinite(vals);
    noEvent  = ~hasEvent;

    if any(hasEvent)
        scatter(P.lon0(hasEvent), P.lat0(hasEvent), 18, vals(hasEvent), 'filled');
        cb = colorbar;
        ylabel(cb, cbLabel);
    end

    if any(noEvent)
        scatter(P.lon0(noEvent), P.lat0(noEvent), 20, 'k', 'x');
    end

    xlabel('Longitude');
    ylabel('Latitude');
    title(strrep(valueName, '_', ' '), 'Interpreter', 'none');

    grid on;
    axis tight;
    box on;

    saveas(fig, outFile);
end

function plotEventHistogram(x, xLabelText, outFile)

    x = x(isfinite(x));

    fig = figure('Color', 'w');

    if isempty(x)
        text(0.5, 0.5, 'No valid values', 'HorizontalAlignment', 'center');
        axis off;
    else
        histogram(x, 40);
        xlabel(xLabelText);
        ylabel('Particle count');
        grid on;
    end

    saveas(fig, outFile);
end

function plotRegionalBoxplot(P, valueName, yLabelText, outFile)

    vals = P.(valueName);
    valid = isfinite(vals);

    fig = figure('Color', 'w');

    if ~any(valid)
        text(0.5, 0.5, 'No valid event times available', ...
            'HorizontalAlignment', 'center');
        axis off;
    else
        boxchart(categorical(string(P.source_region(valid))), vals(valid));
        ylabel(yLabelText);
        xlabel('Source region');
        title(strrep(valueName, '_', ' '), 'Interpreter', 'none');
        grid on;
    end

    saveas(fig, outFile);
end

function plotAnnualMetricByRegion(S, metricName, yLabelText, outFile)

    fig = figure('Color', 'w');
    hold on;

    regions = unique(string(S.group));

    for i = 1:numel(regions)

        r = regions(i);

        idx = string(S.group) == r & isfinite(S.release_year);

        if ~any(idx)
            continue;
        end

        yy = S.release_year(idx);
        val = S.(metricName)(idx);

        [yy, order] = sort(yy);
        val = val(order);

        plot(yy, val, '-o', 'LineWidth', 1.5, 'DisplayName', r);
    end

    xlabel('Release year');
    ylabel(yLabelText);
    title(strrep(metricName, '_', ' '), 'Interpreter', 'none');
    legend('Location', 'best', 'Interpreter', 'none');
    grid on;
    box on;

    saveas(fig, outFile);
end
