%% HKS_cluster_QC_unifiedV7.m
% Batch per-recording cluster QC and spike export. Despite the name this writes
% the complete per-recording preprocessing artifact, not just a mask: downstream
% scripts load this one file instead of redoing NI+AP sync, loadKSdir, depth
% extraction and stimfile parsing.
%
% ---------------------------------------------------------------------
% INPUTS
% ---------------------------------------------------------------------
%   recSpec                    {mouseID, 'all' | {dates}} rows. Folders are found
%                              by globbing catgt* and parsing YYYY-MM-DD. One
%                              recording per day, so the DATE is the unique key --
%                              R01 is not a run index. Each recID must resolve to
%                              exactly one catgt folder; duplicates are an error.
%   <ksdir>                    Kilosort4 output, read via loadKSdir.
%   <ksdir>/cluster_bc_unitType.tsv
%                              Bombcell labels, joined on cluster_id. Columns
%                              cluster_id, bc_unitType, wfDuration_us. A TSV
%                              without wfDuration_us leaves cellType "" and warns.
%   <homedir>/Rec*ProbeLocation*.mat
%                              UniversalProbeFinder anatomy. The 'Rec' prefix is
%                              REQUIRED: a bare *ProbeLocation*.mat matches a
%                              deprecated file and loads the wrong anatomy.
%   <homedir>/*_stimfile_new.mat, <ksdir>/syncSY.mat, <homedir>/*nidq.bin, *ap.meta
%
% ---------------------------------------------------------------------
% OUTPUT  ->  <homedir>/<recName>_cluster_QC.mat   (-v7.3, -struct save)
% ---------------------------------------------------------------------
%   sQC.clusterTable  one row per cluster:
%       clusterID, globalID, nSpikes, depth_um, depthFromSurface, probeDepth_um,
%       ksGroup, bcType, wfDuration_us (us), cellType ("NS"/"BS"/""),
%       firstSpikeLat_ms, peakLat_ms, zetaP_pulse, zetaZ_pulse, meanP_pulse,
%       region, regionAllen, v1Layer (L1..L6b; "" outside V1), recDate, recID,
%       mouseID, lowbandExcess (spk/cycle), responseClass (+1/0/-1),
%       the per-criterion booleans pass_ksLabel, pass_bcLabel, pass_notArtefact,
%       pass_responsive, the informational flags flag_earlyLatency and
%       flag_superficial, and includeDefault
%   sQC.spikes        long-format spike list for the whole recording:
%                     spikeTime (AP clock, rate-corrected), clusterID
%   sQC.stim          pulseOn/Off_AP, flickerOn/Off_AP, flickerFreq (Hz),
%                     dblT0_Correction. Present so downstream scripts never
%                     recompute onsets from structEP.
%   sQC.meta          every threshold and choice used, including dblT0_Correction,
%                     dblSampRateCorrectionImec, ksGoodGroups, bcGoodTypes,
%                     NS_CUTOFF_US, includeGates and regionSource
%   Figures           <recName>_cluster_QC.png / .fig
%
% ---------------------------------------------------------------------
% CURRENT BEHAVIOUR
% ---------------------------------------------------------------------
%   * There is no authoritative include column. includeDefault is a convenience
%     AND of the pass_* flags; inclusion is decided at LOAD time by the
%     downstream toggles, so a gate can be relaxed without re-running QC.
%   * The flat pulse-latency artefact screen is AUDIT-ONLY
%     (cfg.INCLUDE_ARTEFACT_GATE = false). It is rate-confounded, catches
%     essentially no real artefacts, and discarded good units. firstSpikeLat_ms
%     and pass_notArtefact are still computed and stored.
%   * Anatomy comes from a direct cluster_id join against the ProbeLocation file,
%     not from depth bounds. The AP.mat depth-band path is a fallback only, and
%     compares in the SURFACE frame (depthFromSurface).
%   * v1Layer is stored at FINE granularity (L1, L2/3, L4, L5, L6a, L6b) so that
%     grouping -- merging L6a+L6b, dropping L1 -- stays a plot-time decision.
%   * cellType is derived at load time from wfDuration_us at cfg.NS_CUTOFF_US.
%     The width is stored continuous and the label is a derivation. The cutoff is
%     cortex-calibrated, so the label is only interpreted within V1 downstream,
%     and it never enters includeDefault.
%   * Each recording runs inside try/catch, so one failure does not abort the
%     batch; a status table prints at the end. There is no skip-if-exists guard,
%     so every recording in recSpec is recomputed and overwritten on each run.
%   * ZETA resampling is unseeded (250 resamples per cluster), so zetaP_pulse
%     shifts slightly between runs and borderline pass_responsive cells can flip.
%
% ---------------------------------------------------------------------
% RUN ORDER    HKS_bombcell_batchV3 -> THIS -> HKS_pool_resultantsV3
% DEPENDS ON   Acquipix, zetatest, loadKSdir (cortex-lab/spikes), npy-matlab
% CONVENTIONS  fullfile only, no cd; SETUP -> COMPUTE -> PLOTS; -v7.3
% ---------------------------------------------------------------------

clear; close all
set(0, 'DefaultFigureWindowStyle', 'docked')

%% =================== BLOCK 0: USER CONFIG ===================
DRIVE = 'W:';                              % drive letter for the data root and the output tree
addpath(genpath(fullfile(DRIVE, 'Kersten', 'Packages_and_Resources')))

% Recordings to include, one row per mouse. ONE recording per day, so the DATE
% (YYYY-MM-DD) is the unique key; 'R01' is NOT a run index. Value is either:
%   'all'                         -> every catgt* recording for the mouse
%   {'2022-04-12','2022-04-13'}   -> only those recording DATES
recSpec = { ...
    5, 'all'; ...
    6, 'all'; ...
    7, 'all'; ...
    8, 'all'; ...
    10, 'all'};

ksSubfolder = 'kilosort4';

% --- includeDefault composition (the convenience mask; NOT authoritative) ---
cfg.INCLUDE_ARTEFACT_GATE   = false;   % latency screen is audit-only (see header)
cfg.INCLUDE_RESPONSIVE_GATE = true;    % keep ZETA-on-pulse responsiveness in the mask
                                       % (set false for phase-locking-only inclusion;
                                       %  the BC label is the real gate there)

% --- figures ---
cfg.MAKE_FIGURES          = true;      % per-recording 2x2 QC diagnostic figure
cfg.CLOSE_FIGS_AFTER_SAVE = true;      % close after saving (batch: avoid clutter)

% --- (1) Latency artefact screen (AUDIT-ONLY) ---
cfg.LATENCY_THRESHOLD_MS = 20;    % first-spike latency faster than this = suspect
cfg.RESPONSE_WINDOW_S    = 0.05;  % window after pulse to look for the first spike
cfg.SUPERFICIAL_DEPTH_UM = 200;   % secondary (depth) flag only

% --- (2) Responsiveness (ZETA on pulse) ---
cfg.zetaAlpha     = 0.05;         % pass_responsive := zetaP_pulse < zetaAlpha
cfg.intResampNum  = 250;          % ZETA bootstrap iterations
cfg.intLatencyPks = 4;            % ask zetatest for latency estimates too

% --- (3) Curation labels ---
cfg.ksGoodGroups = [2];           % KS4: 0=noise 1=mua 2=good 3=unsorted; [1 2] adds mua
cfg.bcGoodTypes  = {'good'};      % Bombcell types that pass (case-insensitive)
cfg.NS_CUTOFF_US = 400;           % peak-trough width; < = NS, >= = BS. Cortex-calibrated
                                  %     (= Bombcell templateDuration_Ctx_threshold); the label
                                  %     is only interpreted within V1. Descriptor mode
                                  %     'v1celltype' must use the SAME value.

% --- (4) Low-band response classification (driven/suppressed/unaffected) ---
cfg.CLASS_FREQ_LOW   = 2;         % Hz; flicker band used to build the tuning curve
cfg.CLASS_FREQ_HIGH  = 243;       % Hz
cfg.dblAnchorLow     = 4;         % Hz; anchor band low edge
cfg.dblAnchorHigh    = 16;        % Hz; anchor band high edge
cfg.dblAnchorFloor   = 0;         % spk/cycle; 0 = sign only
cfg.dblSpontWindow_s = 3.0;       % pre-trial blank window; matched across QC, pool and raw-spikes

cfg.DRIVE = DRIVE;

%% =================== BLOCK 1: RESOLVE RECORDING FOLDERS ===================
% recID (parsed date+R0N) is the join key. The Bombcell TSV is read from this
% folder's kilosort4 dir, where the Bombcell batch writes it. Inspect the
% printed manifest before letting BLOCK 2 run.

recList = struct('mouseID', {}, 'recName', {}, 'recID', {}, 'recDate', {}, ...
                 'homedir', {}, 'ksdir', {}, 'bcLabelFile', {});

for i = 1:size(recSpec, 1)
    mouseID = recSpec{i, 1};
    want    = recSpec{i, 2};
    ephysDir = fullfile(DRIVE, 'Gulbinaite', 'Resonance_neuropixels', ...
                        'Data_collection', 'Experiment_1', ...
                        sprintf('M%d', mouseID), 'Ephys');

    d = dir(fullfile(ephysDir, 'catgt*'));
    d = d([d.isdir]);
    if isempty(d)
        warning('M%d: no catgt* folders in %s', mouseID, ephysDir); continue;
    end
    names = {d.name};
    nD    = numel(names);

    % parse date + recID (same regexes as the rest of the pipeline)
    recDates = repmat("", nD, 1);
    recIDs   = repmat("", nD, 1);
    for k = 1:nD
        tD = regexp(names{k}, '(\d{4}-\d{2}-\d{2})',     'tokens', 'once');
        tR = regexp(names{k}, '(\d{4}-\d{2}-\d{2}R\d+)', 'tokens', 'once');
        if ~isempty(tD), recDates(k) = string(tD{1}); end
        if ~isempty(tR), recIDs(k) = string(tR{1}); else, recIDs(k) = string(names{k}); end
    end
    hasDate = recDates ~= "";
    if any(~hasDate)
        warning('M%d: %d catgt* folder(s) lack a parseable date - ignored.', ...
                mouseID, sum(~hasDate));
    end

    % date filter ('all' keeps everything; else match recDates against want)
    isAll = (ischar(want) && strcmpi(want,'all')) || ...
            (isstring(want) && isscalar(want) && strcmpi(want,"all"));
    if isAll
        keepMask = hasDate;
    else
        wantDates = string(want);
        keepMask  = hasDate & ismember(recDates, wantDates);
        missing   = setdiff(wantDates, recDates(hasDate));
        if ~isempty(missing)
            warning('M%d: requested dates not found: %s', mouseID, strjoin(missing, ', '));
        end
    end

    % one entry per unique recID. A recID must resolve to exactly one catgt
    % folder; more than one means a stale duplicate is still on disk, which
    % would otherwise be picked silently.
    uID = unique(recIDs(keepMask), 'stable');
    for u = 1:numel(uID)
        sel = find(recIDs == uID(u) & keepMask);
        assert(isscalar(sel), ...
            'M%d %s: %d catgt folders share this recID (%s). Remove the stale one.', ...
            mouseID, uID(u), numel(sel), strjoin(names(sel), ', '));
        qcIdx = sel;

        home  = fullfile(ephysDir, names{qcIdx});
        ksdir = fullfile(home, ksSubfolder);
        if ~exist(ksdir, 'dir')
            warning('M%d %s: no %s under %s - skipping.', ...
                    mouseID, uID(u), ksSubfolder, home); continue;
        end

        % Bombcell TSV is written into this folder's kilosort4 dir
        bcFile = fullfile(ksdir, 'cluster_bc_unitType.tsv');

        e = numel(recList) + 1;
        recList(e).mouseID     = mouseID;
        recList(e).recName     = names{qcIdx};
        recList(e).recID       = char(uID(u));
        recList(e).recDate     = char(recDates(sel(1)));
        recList(e).homedir     = home;
        recList(e).ksdir       = ksdir;
        recList(e).bcLabelFile = bcFile;
    end
end

% chronological order across the whole set (by date, then mouse)
if ~isempty(recList)
    [~, ord] = sortrows([string({recList.recDate})', [recList.mouseID]']);
    recList  = recList(ord);
end

nRec = numel(recList);
fprintf('\nResolved %d recordings (chronological):\n', nRec);
for e = 1:nRec
    bcTag = 'BC-tsv ok'; if ~exist(recList(e).bcLabelFile, 'file'), bcTag = 'BC-tsv MISSING'; end
    fprintf('  M%-2d  %-14s  %-40s  [%s]\n', ...
            recList(e).mouseID, recList(e).recID, recList(e).recName, bcTag);
end
if nRec == 0, error('No recordings resolved - check DRIVE and recSpec.'); end

%% =================== BLOCK 2: RUN QC PER RECORDING ===================
summary = struct('recID', {}, 'mouseID', {}, 'status', {}, 'nClust', {}, ...
                 'nKS', {}, 'nBC', {}, 'nNotArt', {}, 'nResp', {}, ...
                 'nInclude', {}, 'regionSource', {}, 'msg', {});

for e = 1:nRec
    r = recList(e);
    fprintf('\n====================================\n');
    fprintf('[%d/%d] M%d  %s  (%s)\n', e, nRec, r.mouseID, r.recID, r.recName);
    fprintf('====================================\n');
    try
        s = run_qc_one(r, cfg);
        s.status = 'ok'; s.msg = '';
    catch ME
        s = struct('recID', r.recID, 'mouseID', r.mouseID, 'status', 'FAILED', ...
                   'nClust', NaN, 'nKS', NaN, 'nBC', NaN, 'nNotArt', NaN, ...
                   'nResp', NaN, 'nInclude', NaN, 'regionSource', 'n/a', ...
                   'msg', ME.message);
        fprintf(2, '  FAILED: %s\n', ME.message);
    end
    summary(e) = s; %#ok<SAGROW>
end

%% =================== BLOCK 3: BATCH SUMMARY ===================
fprintf('\n================= BATCH SUMMARY =================\n');
T = struct2table(summary);
disp(T(:, {'mouseID','recID','status','nClust','nKS','nBC', ...
           'nNotArt','nResp','nInclude','regionSource'}));

nFail = sum(strcmp({summary.status}, 'FAILED'));
if nFail > 0
    fprintf(2, '\n%d recording(s) FAILED:\n', nFail);
    for e = 1:nRec
        if strcmp(summary(e).status, 'FAILED')
            fprintf(2, '  M%d %s: %s\n', summary(e).mouseID, summary(e).recID, summary(e).msg);
        end
    end
else
    fprintf('All %d recordings completed.\n', nRec);
end

%% ==========================================================================
%% HOW TO LOAD A QC FILE IN A PLOTTING / ANALYSIS SCRIPT
%% ==========================================================================
% sQC = load(fullfile(homedir, sprintf('%s_cluster_QC.mat', recName)));
% C   = sQC.clusterTable;
%
% % --- Choose which QC gates to enforce (toggle freely) ---
% % These are the settings the pooled analysis was built with. Changing them
% % changes the cluster set and will not reproduce the reported numbers.
% useKS         = false;   % Kilosort 'group' proved unreliable against Bombcell;
%                          % pass_bcLabel is the authoritative label gate
% useBC         = true;    % require Bombcell 'good'
% useResponsive = true;    % pulse-ZETA responders only. Set false for a
%                          % phase-locking analysis: a rate-unaffected cell can
%                          % still phase-lock, and excluding it begs the question
% useArtefact   = false;   % flat pulse-latency screen is an audit column, not a
%                          % gate (rate-confounded; ~0 real artefacts)
%
% incl = true(height(C), 1);
% if useKS,         incl = incl & C.pass_ksLabel;     end
% if useBC,         incl = incl & C.pass_bcLabel;     end
% if useResponsive, incl = incl & C.pass_responsive;  end
% if useArtefact,   incl = incl & C.pass_notArtefact; end
%
% keepClusterIDs = C.clusterID(incl);
%
% % --- Pull spikes for kept clusters from the long-format list ---
% S   = sQC.spikes;
% sel = ismember(S.clusterID, keepClusterIDs);
% keptTimes = S.spikeTime(sel);  keptClusters = S.clusterID(sel);

%% =================== LOCAL FUNCTIONS ===================

function s = run_qc_one(r, cfg)
% One recording's full QC pass. Writes <recName>_cluster_QC.mat into r.homedir
% and returns a one-row summary struct.
mouseID = r.mouseID;  recName = r.recName;  recID = r.recID;  recDate = r.recDate;
homedir = r.homedir;  ksdir   = r.ksdir;    bcLabelFile = r.bcLabelFile;

s = struct('recID', recID, 'mouseID', mouseID, 'status', '', 'nClust', 0, ...
           'nKS', 0, 'nBC', 0, 'nNotArt', 0, 'nResp', 0, 'nInclude', 0, ...
           'regionSource', 'none', 'msg', '');

strFigDir = fullfile(cfg.DRIVE, 'Kersten', 'Figures', sprintf('M%d', mouseID));
if cfg.MAKE_FIGURES && ~isempty(strFigDir) && ~exist(strFigDir, 'dir'), mkdir(strFigDir); end

%% ---- NI timing (calibrated rate + T0) ----
ni_files = dir(fullfile(homedir, '*nidq.bin'));
if isempty(ni_files), error('No *nidq.bin in %s', homedir); end
sMetaNI   = DP_ReadMeta(ni_files(1).name, homedir);
matDataNI = -DP_ReadBin(-inf, inf, sMetaNI, ni_files(1).name, homedir);

intCheckSyncPulseCh = PP_FindSyncPulseCh(matDataNI);
boolVecSyncPulsesNI = DP_GetUpDown(matDataNI(intCheckSyncPulseCh, :));
dblSampRateNI       = PP_GetPulseIntervalFromBinVec(boolVecSyncPulsesNI);
intFirstSampleNI    = str2double(sMetaNI.firstSample);
dblT0_NI_new        = intFirstSampleNI / dblSampRateNI;
clear matDataNI boolVecSyncPulsesNI

%% ---- AP timing (calibrated rate + T0) from syncSY.mat ----
ap_files = dir(fullfile(homedir, '*ap.meta'));
if isempty(ap_files), error('No *ap.meta in %s', homedir); end
sMetaAP = DP_ReadMeta(ap_files(1).name, homedir);

syncFile = fullfile(ksdir, 'syncSY.mat');
if ~exist(syncFile, 'file'), syncFile = fullfile(homedir, 'syncSY.mat'); end
if ~exist(syncFile, 'file')
    error('syncSY.mat not found in %s or %s. Run Acquipix preprocessing first.', ksdir, homedir);
end
sSync = load(syncFile);

vecChangePulsesImec = diff(sSync.syncSY);
vecSyncPulseOnImec  = find(vecChangePulsesImec == 1) + 1;
vecDiffPulsesImec   = sort(diff(vecSyncPulseOnImec));
intNumPeriods       = numel(vecDiffPulsesImec);
intLo               = ceil(intNumPeriods/10);
intHi               = floor(intNumPeriods*9/10);
dblSampRateImAp     = mean(vecDiffPulsesImec(intLo:intHi));
intFirstSampleImAp  = str2double(sMetaAP.firstSample);
dblT0_ImAp_new      = intFirstSampleImAp / dblSampRateImAp;
dblT0_Correction    = dblT0_ImAp_new - dblT0_NI_new;
fprintf('T0 correction: %.3f ms\n', dblT0_Correction * 1000);

%% ---- Load Kilosort spikes + apply rate correction ----
sSpikes = loadKSdir(ksdir);
dblImecSampRateReported   = sSpikes.sample_rate;
dblSampRateCorrectionImec = dblImecSampRateReported / dblSampRateImAp;

cids   = sSpikes.cids(:);
nClust = numel(cids);
cellSpikes = cell(nClust, 1);
for ii = 1:nClust
    cellSpikes{ii} = sSpikes.st(sSpikes.clu == cids(ii)) * dblSampRateCorrectionImec;
end
fprintf('Loaded %d clusters.\n', nClust);

%% ---- Per-cluster depth from Kilosort templates ----
[~, spikeDepths] = templatePositionsAmplitudes(sSpikes.temps, sSpikes.winv, ...
    sSpikes.ycoords, sSpikes.spikeTemplates, sSpikes.tempScalingAmps);
clusterDepths = nan(nClust, 1);
for ii = 1:nClust
    clusterDepths(ii) = nanmedian(spikeDepths(sSpikes.clu == cids(ii)));
end
probeMaxY        = max(sSpikes.ycoords);       % top channel = surface
depthFromSurface = probeMaxY - clusterDepths;  % 0 = surface, large = deep

%% ---- Anatomical region per cluster (cluster_id join, ProbeLocation) ----
% Area labels come from the UPF ProbeLocation file by a DIRECT cluster_id join
% (0-indexed, matches sSpikes.cids), NOT AP.mat depth bounds. The AP-band
% fallback below compares in the SURFACE frame (depthFromSurface).
% File glob MUST keep the 'Rec' prefix (a bare *ProbeLocation*.mat matches a
% deprecated file and would load the wrong anatomy).
areaNames     = {'V1',                   'LGN',                           'CA1',       'CA3'};
areaQueries   = {'Primary visual area', 'Dorsal part of the lateral ',   'Field CA1', 'Field CA3'};
region        = repmat("", nClust, 1);
regionAllen   = repmat("", nClust, 1);
v1Layer       = repmat("", nClust, 1);    % V1 laminar label; "" outside V1 (ProbeLocation only)
probeDepth_um = nan(nClust, 1);
areaBounds    = nan(numel(areaNames), 2);
regionSource  = 'none';

plFile = dir(fullfile(homedir, 'Rec*ProbeLocation*.mat'));
if ~isempty(plFile)
    regionSource = 'ProbeLocation';
    sPL = load(fullfile(homedir, plFile(1).name));
    if ~isfield(sPL, 'sProbeCoords') || ~isfield(sPL.sProbeCoords, 'sProbeAdjusted')
        error('ProbeLocation file %s lacks sProbeCoords.sProbeAdjusted.', plFile(1).name);
    end
    sPA = sPL.sProbeCoords.sProbeAdjusted;
    for f = {'cluster_id', 'probe_area_full_per_cluster'}
        if ~isfield(sPA, f{1})
            error('ProbeLocation sProbeAdjusted is missing required field "%s".', f{1});
        end
    end
    plClusterID = double(sPA.cluster_id(:));
    plAreaFull  = string(sPA.probe_area_full_per_cluster(:));
    [tf, loc]   = ismember(cids, plClusterID);
    regionAllen(tf) = plAreaFull(loc(tf));
    if isfield(sPA, 'depth_per_cluster')
        plDepth = double(sPA.depth_per_cluster(:));
        probeDepth_um(tf) = plDepth(loc(tf));
    end
    for a = 1:numel(areaNames)
        region(contains(regionAllen, areaQueries{a})) = string(areaNames{a});
    end
    % V1 laminar label (Schneider-style), parsed from the verbatim Allen name
    % ("Primary visual area, layer 2/3" -> "L2/3"). Coarse `region` is untouched;
    % fine layers (L1, L2/3, L4, L5, L6a, L6b) are kept as-is so the grouping
    % (e.g. 6a+6b -> 6, drop/keep L1) can be chosen at plot time. Populated for
    % V1 units only; "" elsewhere. (AP-legacy fallback leaves this "".)
    for ii = find(region == "V1")'
        t = regexp(regionAllen(ii), 'layer\s+([\w/]+)', 'tokens', 'once', 'ignorecase');
        if ~isempty(t), v1Layer(ii) = "L" + string(t{1}); end
    end
    if any(region == "V1")
        uL = unique(v1Layer(region == "V1"));
        fprintf('V1 layers:');
        for k = 1:numel(uL)
            lab = uL(k); if lab == "", lab = "(none)"; end
            fprintf(' %s=%d', lab, sum(region == "V1" & v1Layer == uL(k)));
        end
        fprintf('\n');
    end
    matchFrac = mean(tf);
    fprintf('ProbeLocation join: %d / %d clusters matched cluster_id (%.1f%%) from %s\n', ...
            sum(tf), nClust, 100*matchFrac, plFile(1).name);
    if matchFrac < 0.95
        warning(['Only %.1f%% of KS4 cids matched ProbeLocation cluster_id - usually a ' ...
                 'sort mismatch. Verify sort provenance before trusting the labels.'], 100*matchFrac);
    end
    fprintf('Region (coarse): V1=%d LGN=%d CA1=%d CA3=%d  (no coarse area: %d of %d matched)\n', ...
            sum(region=="V1"), sum(region=="LGN"), sum(region=="CA1"), ...
            sum(region=="CA3"), sum(tf & region==""), sum(tf));
else
    warning(['No Rec*ProbeLocation*.mat in %s. Falling back to the legacy AP.mat ' ...
             'depth-band path (less exact at borders).'], homedir);
    APfile = dir(fullfile(homedir, '*AP.mat'));
    if isempty(APfile)
        warning('No *AP.mat either; region = "" for all clusters.');
    else
        sAPdat = load(fullfile(homedir, APfile(1).name));
        if ~isfield(sAPdat, 'sAP')
            warning('AP.mat has no sAP variable; region left "".');
        else
            regionSource = 'AP-legacy';
            sAPAreas    = {sAPdat.sAP.sCluster.Area};
            sAPDepthAll = [sAPdat.sAP.sCluster.Depth];      % SURFACE-referenced
            for a = 1:numel(areaNames)
                inArea = contains(sAPAreas, areaQueries{a});
                dd = sAPDepthAll(inArea); dd = dd(~isnan(dd));
                if ~isempty(dd), areaBounds(a,:) = [min(dd), max(dd)]; end
            end
            nOverlap = 0;
            for ii = 1:nClust
                hits = find(depthFromSurface(ii) >= areaBounds(:,1) & ...
                            depthFromSurface(ii) <= areaBounds(:,2));
                if isempty(hits), continue; end
                if numel(hits) > 1
                    nOverlap = nOverlap + 1;
                    [~, k] = min(abs(mean(areaBounds(hits,:),2) - depthFromSurface(ii)));
                    hits = hits(k);
                end
                region(ii) = string(areaNames{hits});
            end
            if nOverlap > 0
                warning('%d clusters in overlapping bands; assigned to nearest-centre area.', nOverlap);
            end
            fprintf('Region (AP-legacy): V1=%d LGN=%d CA1=%d CA3=%d unassigned=%d\n', ...
                    sum(region=="V1"), sum(region=="LGN"), sum(region=="CA1"), ...
                    sum(region=="CA3"), sum(region==""));
        end
    end
end

%% ---- Stimulus structure: PULSE onsets on the AP clock ----
stim_files = dir(fullfile(homedir, '*_stimfile_new.mat'));
if isempty(stim_files), error('No *_stimfile_new.mat found in %s', homedir); end
load(fullfile(homedir, stim_files(1).name));   % structEP

pulse_idx   = find(contains(structEP.vecTrialStimTypes, 'Pulse'));
pulseOn_AP  = structEP.ActOnNI_uncorr(pulse_idx)  - dblT0_Correction;
pulseOff_AP = structEP.ActOffNI_uncorr(pulse_idx) - dblT0_Correction;
fprintf('Pulse trials: %d\n', numel(pulseOn_AP));
if isempty(pulseOn_AP)
    error('No Pulse trials found - the latency + ZETA screens both need pulse onsets.');
end

%% ---- Bombcell unit-type labels ----
bcType   = repmat({''}, nClust, 1);
wfDur_us = nan(nClust, 1);                 % peak-trough width, us; NaN if the TSV lacks the column
if exist(bcLabelFile, 'file')
    bcTbl = readtable(bcLabelFile, 'FileType', 'text', 'Delimiter', '\t');
    [tf, loc] = ismember(cids, bcTbl.cluster_id);
    bcType(tf) = bcTbl.bc_unitType(loc(tf));
    fprintf('Bombcell labels: %d / %d clusters matched from %s\n', sum(tf), nClust, bcLabelFile);
    if sum(tf) < nClust
        warning('%d clusters have no Bombcell label (pass_bcLabel=false for those).', nClust - sum(tf));
    end
    % waveform width rides the SAME join as bc_unitType (no separate ismember)
    if ismember('wfDuration_us', bcTbl.Properties.VariableNames)
        wfDur_us(tf) = bcTbl.wfDuration_us(loc(tf));
    else
        warning(['TSV %s has no wfDuration_us column. ' ...
                 'cellType will be "" for this recording -- rerun HKS_bombcell_batchV3 ' ...
                 'to populate it.'], bcLabelFile);
    end
else
    warning(['Bombcell label file not found at %s. pass_bcLabel = false for ALL clusters ' ...
             '- leave the BC toggle OFF at load time.'], bcLabelFile);
end

% NS/BS from width. NaN width -> "" (unlabelled -> excluded from cell-type
% grouping downstream). Cortex-calibrated cutoff; only interpreted within V1.
cellType = repmat("", nClust, 1);
cellType(wfDur_us <  cfg.NS_CUTOFF_US) = "NS";
cellType(wfDur_us >= cfg.NS_CUTOFF_US) = "BS";

%% ---- Per-cluster QC metrics: first-spike latency + ZETA on pulse ----
matStimOnOff = [pulseOn_AP(:), pulseOff_AP(:)];
if numel(pulseOn_AP) > 1
    dblUseMaxDur = nanmin(diff(matStimOnOff(:,1)));
else
    dblUseMaxDur = nanmedian(pulseOff_AP - pulseOn_AP);
end
if ~isfinite(dblUseMaxDur) || dblUseMaxDur <= 0, dblUseMaxDur = cfg.RESPONSE_WINDOW_S; end
intPlot = 0;

firstSpikeLat_ms = nan(nClust, 1);
peakLat_ms       = nan(nClust, 1);
zetaP_pulse      = nan(nClust, 1);
zetaZ_pulse      = nan(nClust, 1);
meanP_pulse      = nan(nClust, 1);

fprintf('Running latency + ZETA screen on %d clusters (pulse)...\n', nClust);
tic
for ii = 1:nClust
    spikes = cellSpikes{ii};
    if numel(spikes) < 10, continue; end

    % ZETA against pulse onsets
    try
        [zP, sZETA] = zetatest(spikes, matStimOnOff, dblUseMaxDur, ...
                               cfg.intResampNum, intPlot, cfg.intLatencyPks, [], [], false);
        zetaP_pulse(ii) = zP;
        if ~isempty(sZETA) && isstruct(sZETA)
            if isfield(sZETA, 'dblZETA'),  zetaZ_pulse(ii) = sZETA.dblZETA; end
            if isfield(sZETA, 'dblMeanP'), meanP_pulse(ii) = sZETA.dblMeanP; end
            if isfield(sZETA, 'vecLatencies') && ~isempty(sZETA.vecLatencies)
                peakLat_ms(ii) = sZETA.vecLatencies(1) * 1000;
            end
        end
    catch
    end

    % first-spike latency per pulse trial, then median
    lats = nan(numel(pulseOn_AP), 1);
    for t = 1:numel(pulseOn_AP)
        post = spikes(spikes > pulseOn_AP(t) & spikes < pulseOn_AP(t) + cfg.RESPONSE_WINDOW_S);
        if ~isempty(post), lats(t) = (post(1) - pulseOn_AP(t)) * 1000; end
    end
    firstSpikeLat_ms(ii) = nanmedian(lats);
end
fprintf('  done in %.1f s.\n', toc);

%% ---- Flicker trials on the AP clock (for classification + sQC.stim) ----
isFlicker     = ~ismember(structEP.vecTrialStimTypes, {'Sparse','Checker','Pulse','WN'});
flickerLabels = structEP.vecTrialStimTypes(isFlicker);
flickerOn_AP  = structEP.ActOnNI_uncorr(isFlicker)  - dblT0_Correction;
flickerOff_AP = structEP.ActOffNI_uncorr(isFlicker) - dblT0_Correction;
flickerFreq   = cellfun(@(s) sscanf(s, '%f Hz'), flickerLabels);
isInBand      = flickerFreq >= cfg.CLASS_FREQ_LOW & flickerFreq <= cfg.CLASS_FREQ_HIGH & ~isnan(flickerFreq);
flickerOn_AP  = flickerOn_AP(isInBand);  flickerOff_AP = flickerOff_AP(isInBand);
flickerFreq   = flickerFreq(isInBand);   flickerFreq   = flickerFreq(:);
nFlickTrials  = numel(flickerOn_AP);
flickerDur    = flickerOff_AP(:) - flickerOn_AP(:);
flickerNCyc   = flickerDur .* flickerFreq;
vecUniqFreq   = unique(flickerFreq);
fprintf('Flicker trials in [%g,%g] Hz: %d\n', cfg.CLASS_FREQ_LOW, cfg.CLASS_FREQ_HIGH, nFlickTrials);

%% ---- Build per-criterion pass flags ----
% (1) Not a latency artefact (AUDIT-ONLY: not folded into includeDefault).
flag_earlyLatency = firstSpikeLat_ms < cfg.LATENCY_THRESHOLD_MS;   % NaN -> false
pass_notArtefact  = ~flag_earlyLatency;
% Secondary depth flag (informational).
flag_superficial  = depthFromSurface < cfg.SUPERFICIAL_DEPTH_UM;
% (2) Responsive to pulse by ZETA.
pass_responsive   = zetaP_pulse < cfg.zetaAlpha;                   % NaN -> false
% (3a) Kilosort label.
ksGroup      = sSpikes.cgs(:);
pass_ksLabel = ismember(ksGroup, cfg.ksGoodGroups);
% (3b) Bombcell label (case-insensitive; missing label fails).
pass_bcLabel = false(nClust, 1);
for ii = 1:nClust
    if ~isempty(bcType{ii}), pass_bcLabel(ii) = any(strcmpi(bcType{ii}, cfg.bcGoodTypes)); end
end

%% ---- Low-band response classification ----
lowbandExcess = nan(nClust, 1);
responseClass = zeros(nClust, 1);
if nFlickTrials == 0
    warning('No flicker trials; responseClass = 0 for all.');
else
    matCount = zeros(nClust, nFlickTrials);
    for ii = 1:nClust
        spk = cellSpikes{ii};
        if isempty(spk), continue; end
        for t = 1:nFlickTrials
            matCount(ii,t) = sum(spk >= flickerOn_AP(t) & spk < flickerOff_AP(t));
        end
    end
    vecSpontRate       = computeSpontRate(cellSpikes, flickerOn_AP, cfg.dblSpontWindow_s);
    matSpkPerCycleCorr = (matCount - vecSpontRate .* flickerDur(:)') ./ flickerNCyc(:)';
    matMu3 = nan(nClust, numel(vecUniqFreq));
    for fi = 1:numel(vecUniqFreq)
        sel = flickerFreq == vecUniqFreq(fi);
        if any(sel), matMu3(:,fi) = mean(matSpkPerCycleCorr(:,sel), 2, 'omitnan'); end
    end
    inAnchor      = vecUniqFreq >= cfg.dblAnchorLow & vecUniqFreq <= cfg.dblAnchorHigh;
    lowbandExcess = mean(matMu3(:, inAnchor), 2, 'omitnan');
    responseClass = double(lowbandExcess >  cfg.dblAnchorFloor) - ...
                    double(lowbandExcess < -cfg.dblAnchorFloor);
    fprintf('Low-band class (anchor %.3g-%.3g Hz, floor %.3g): driven=%d suppr=%d unaff=%d\n', ...
            cfg.dblAnchorLow, cfg.dblAnchorHigh, cfg.dblAnchorFloor, ...
            sum(responseClass>0), sum(responseClass<0), sum(responseClass==0));
end

%% ---- Assemble per-cluster table ----
clusterTable = table();
clusterTable.clusterID         = cids;
clusterTable.globalID          = compose("%s__clu%d", string(recName), cids);
clusterTable.nSpikes           = cellfun(@numel, cellSpikes);
clusterTable.depth_um          = clusterDepths;
clusterTable.depthFromSurface  = depthFromSurface;
clusterTable.ksGroup           = ksGroup;
clusterTable.bcType            = string(bcType);
clusterTable.wfDuration_us     = wfDur_us;            % peak-trough width, us (continuous)
clusterTable.firstSpikeLat_ms  = firstSpikeLat_ms;
clusterTable.peakLat_ms        = peakLat_ms;
clusterTable.zetaP_pulse       = zetaP_pulse;
clusterTable.zetaZ_pulse       = zetaZ_pulse;
clusterTable.meanP_pulse       = meanP_pulse;
% --- per-criterion pass flags (toggleable QC gates) ---
clusterTable.pass_ksLabel      = pass_ksLabel;
clusterTable.pass_bcLabel      = pass_bcLabel;
clusterTable.pass_notArtefact  = pass_notArtefact;    % audit column only
clusterTable.pass_responsive   = pass_responsive;
% --- secondary / informational flags ---
clusterTable.flag_earlyLatency = flag_earlyLatency;
clusterTable.flag_superficial  = flag_superficial;
clusterTable.mouseID       = repmat(mouseID, nClust, 1);
clusterTable.recID         = repmat(string(recID),   nClust, 1);
clusterTable.recDate       = repmat(string(recDate), nClust, 1);
clusterTable.region        = region;
clusterTable.regionAllen   = regionAllen;
clusterTable.v1Layer       = v1Layer;          % V1 laminar label (L1..L6b); "" outside V1
clusterTable.cellType      = cellType;         % NS/BS from wfDuration_us; "" if unlabelled
clusterTable.probeDepth_um = probeDepth_um;
clusterTable.lowbandExcess = lowbandExcess;
clusterTable.responseClass = responseClass;

% Convenience "default include" mask (NOT authoritative). The flat pulse-latency
% artefact gate is excluded unless cfg.INCLUDE_ARTEFACT_GATE is set.
includeDefault = pass_ksLabel & pass_bcLabel;
gateList = {'ksLabel','bcLabel'};
if cfg.INCLUDE_RESPONSIVE_GATE
    includeDefault = includeDefault & pass_responsive; gateList{end+1} = 'responsive';
end
if cfg.INCLUDE_ARTEFACT_GATE
    includeDefault = includeDefault & pass_notArtefact; gateList{end+1} = 'notArtefact';
end
clusterTable.includeDefault = includeDefault;

fprintf('\nPer-criterion pass counts (of %d clusters):\n', nClust);
fprintf('  pass_ksLabel     : %d\n', sum(pass_ksLabel));
fprintf('  pass_bcLabel     : %d\n', sum(pass_bcLabel));
fprintf('  pass_notArtefact : %d  (audit only)\n', sum(pass_notArtefact));
fprintf('  pass_responsive  : %d\n', sum(pass_responsive));
fprintf('  includeDefault [%s] : %d\n', strjoin(gateList, ' & '), sum(includeDefault));

%% ---- Assemble LONG-format spike list ----
vecAllTimes     = vertcat(cellSpikes{:});
vecClustOfSpike = zeros(numel(vecAllTimes), 1);
iStart = 1;
for ii = 1:nClust
    n = numel(cellSpikes{ii});
    if n == 0, continue; end
    vecClustOfSpike(iStart:iStart+n-1) = cids(ii);
    iStart = iStart + n;
end
[vecAllTimes, ord] = sort(vecAllTimes);
vecClustOfSpike    = vecClustOfSpike(ord);

spikes = struct();
spikes.spikeTime = vecAllTimes;
spikes.clusterID = vecClustOfSpike;
spikes.recName   = string(recName);
spikes.mouseID   = mouseID;
fprintf('Long-format spike list: %d spikes across %d clusters.\n', numel(vecAllTimes), nClust);

%% ---- Recording-level metadata ----
meta = struct();
meta.recName                   = recName;
meta.mouseID                   = mouseID;
meta.ksdir                     = ksdir;
meta.homedir                   = homedir;
meta.bcLabelFile               = bcLabelFile;
meta.dblT0_Correction          = dblT0_Correction;
meta.dblSampRateCorrectionImec = dblSampRateCorrectionImec;
meta.dblSampRateImAp           = dblSampRateImAp;
meta.nPulseTrials              = numel(pulseOn_AP);
meta.LATENCY_THRESHOLD_MS      = cfg.LATENCY_THRESHOLD_MS;
meta.RESPONSE_WINDOW_S         = cfg.RESPONSE_WINDOW_S;
meta.SUPERFICIAL_DEPTH_UM      = cfg.SUPERFICIAL_DEPTH_UM;
meta.zetaAlpha                 = cfg.zetaAlpha;
meta.intResampNum              = cfg.intResampNum;
meta.ksGoodGroups              = cfg.ksGoodGroups;
meta.bcGoodTypes               = {cfg.bcGoodTypes};
meta.NS_CUTOFF_US              = cfg.NS_CUTOFF_US;    % NS/BS peak-trough cutoff (us)
meta.INCLUDE_ARTEFACT_GATE     = cfg.INCLUDE_ARTEFACT_GATE;
meta.INCLUDE_RESPONSIVE_GATE   = cfg.INCLUDE_RESPONSIVE_GATE;
meta.includeGates              = gateList;
meta.savedOn                   = datestr(now, 'yyyy-mm-dd HH:MM:SS');
meta.areaNames = areaNames; meta.areaQueries = areaQueries; meta.areaBounds = areaBounds;
meta.regionSource = regionSource;
meta.dblAnchorLow = cfg.dblAnchorLow; meta.dblAnchorHigh = cfg.dblAnchorHigh;
meta.dblAnchorFloor = cfg.dblAnchorFloor; meta.dblSpontWindow_s = cfg.dblSpontWindow_s;
meta.recDate = recDate; meta.recID = recID;

stim = struct();
stim.pulseOn_AP   = pulseOn_AP(:);    stim.pulseOff_AP   = pulseOff_AP(:);
stim.flickerOn_AP = flickerOn_AP(:);  stim.flickerOff_AP = flickerOff_AP(:);
stim.flickerFreq  = flickerFreq(:);   stim.dblT0_Correction = dblT0_Correction;
stim.note = 'AP clock = structEP.ActOnNI_uncorr - dblT0_Correction.';

%% ---- Save unified QC file ----
sQC = struct();
sQC.clusterTable = clusterTable;
sQC.spikes       = spikes;
sQC.stim         = stim;      % pool_resultants reads these instead of reparsing structEP
sQC.meta         = meta;

outFile = fullfile(homedir, sprintf('%s_cluster_QC.mat', recName));
save(outFile, '-struct', 'sQC', '-v7.3');
fprintf('Saved unified QC to %s\n', outFile);

%% ---- Diagnostic figure (optional) ----
if cfg.MAKE_FIGURES
    fh = figure('Name', [recName ' - QC'], 'Color', 'w');
    subplot(2,2,1)
    histogram(firstSpikeLat_ms, 0:1:50); hold on
    plot([cfg.LATENCY_THRESHOLD_MS cfg.LATENCY_THRESHOLD_MS], ylim, 'r--', 'LineWidth', 1.5)
    xlabel('First-spike latency (ms)'); ylabel('# units')
    title(sprintf('Pulse latency (audit cutoff = %d ms)', cfg.LATENCY_THRESHOLD_MS)); box off

    subplot(2,2,2)
    scatter(depthFromSurface(~pass_responsive), firstSpikeLat_ms(~pass_responsive), 20, [.7 .7 .7], 'filled'); hold on
    scatter(depthFromSurface(pass_responsive),  firstSpikeLat_ms(pass_responsive),  20, 'r', 'filled')
    scatter(depthFromSurface(~pass_notArtefact), firstSpikeLat_ms(~pass_notArtefact), 40, 'k', 'LineWidth', 1.2)
    xline(cfg.SUPERFICIAL_DEPTH_UM, 'b--'); yline(cfg.LATENCY_THRESHOLD_MS, 'r--')
    xlabel('Depth from surface (\mum)'); ylabel('First-spike latency (ms)')
    title('Latency vs depth (red=responsive, black ring=audit-artefact)'); box off

    subplot(2,2,3)
    counts = [sum(pass_ksLabel), sum(pass_bcLabel), sum(pass_notArtefact), ...
              sum(pass_responsive), sum(includeDefault)];
    bar(counts); ylim([0 nClust])
    set(gca, 'XTickLabel', {'KS','BC','notArt*','resp','INCL'})
    ylabel('# clusters passing'); title('Per-criterion pass counts (*audit)'); box off

    subplot(2,2,4)
    scatter(firstSpikeLat_ms, zetaZ_pulse, 20, 'k', 'filled'); hold on
    scatter(firstSpikeLat_ms(~pass_notArtefact), zetaZ_pulse(~pass_notArtefact), 30, 'r', 'filled')
    xline(cfg.LATENCY_THRESHOLD_MS, 'r--')
    xlabel('First-spike latency (ms)'); ylabel('ZETA Z (pulse)')
    title('Responsiveness vs latency (red = audit-artefact)'); box off

    if ~isempty(strFigDir)
        savefig(fh, fullfile(strFigDir, sprintf('%s_cluster_QC.fig', recName)));
        saveas( fh, fullfile(strFigDir, sprintf('%s_cluster_QC.png', recName)));
    end
    if cfg.CLOSE_FIGS_AFTER_SAVE, close(fh); end
end

%% ---- summary row ----
s.status = 'ok';  s.nClust = nClust;
s.nKS = sum(pass_ksLabel);  s.nBC = sum(pass_bcLabel);
s.nNotArt = sum(pass_notArtefact);  s.nResp = sum(pass_responsive);
s.nInclude = sum(includeDefault);  s.regionSource = regionSource;  s.msg = '';
fprintf('Done.\n');
end

% ---- per-cluster spontaneous rate (spk/s) --------------------------
% From the pre-trial blank window [onset - dblWindow_s, onset), pooled across
% all supplied onsets.
% Shared helper: duplicated verbatim. Edit every copy or none (see README).
function vecSpontRate = computeSpontRate(cellSpikes, vecOnsets, dblWindow_s)
    intClustNum = numel(cellSpikes);
    intNTrials  = numel(vecOnsets);
    vecSpontRate = zeros(intClustNum, 1);
    vecStart = vecOnsets(:) - dblWindow_s;
    vecEnd   = vecOnsets(:);
    for intClust = 1:intClustNum
        spk = cellSpikes{intClust};
        if isempty(spk), continue; end
        intCount = 0;
        for intTrial = 1:intNTrials
            intCount = intCount + sum(spk >= vecStart(intTrial) & spk < vecEnd(intTrial));
        end
        vecSpontRate(intClust) = intCount / (intNTrials * dblWindow_s);
    end
end
