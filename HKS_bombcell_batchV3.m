%% HKS_bombcell_batchV3.m
% Stage 0: runs Bombcell per recording, pools every unit into one table, applies
% a custom RPV promotion, and writes the per-recording label file the rest of the
% pipeline reads. Also carries the QC-supplement diagnostics.
%
% ---------------------------------------------------------------------
% INPUTS
% ---------------------------------------------------------------------
%   recSpec              {mouseID, 'all' | {dates}} rows. Folders are found by
%                        globbing catgt* and parsing YYYY-MM-DD. One recording per
%                        day, so the DATE is the unique key -- R01 is not a run
%                        index. The date filter is honoured.
%   <homedir>/*ap.bin    Raw AP file. A folder without one is skipped with a warning.
%   <homedir>/*.ap.meta  Read for the Bombcell param defaults.
%   <ksdir>              Kilosort4 output, read via bc.load.loadEphysData.
%
% ---------------------------------------------------------------------
% OUTPUT  ->  <ksdir>/cluster_bc_unitType.tsv          (one per recording)
% ---------------------------------------------------------------------
%   cluster_id     Phy cluster id, the join key used downstream
%   bc_unitType    FINAL label, after RPV promotion. This file is the ONLY thing
%                  this script sends downstream; everything else is diagnostic.
%
% OUTPUT  ->  <DRIVE>\Kersten\QC\bombcell_pooled\bombcell_pooled_QC.mat  (-v7.3)
%   pooled          every unit, keyed by globalID = "<recName>__clu<phyID>";
%                   carries bcLabel (final), bcLabelRaw (pre-promotion),
%                   rpvPromoted, the 9 metricSpec metrics, excl_label,
%                   incl_pctl and keep
%   paramStore      Bombcell's own param struct per recording -- the thresholds
%                   that actually produce the GOOD label, for the Methods table
%   metricPassBC    per-metric pass/fail against BOMBCELL's thresholds (gating)
%   metricPassPctl  per-metric pass/fail against the percentile cutoffs
%                   (descriptive; see CURRENT BEHAVIOUR)
%   passMatBC, passMatPctl
%                   the same as raw logicals, so any other cross-tabulation can
%                   be made later without re-running the Bombcell loop
%   reconTab        label-reconstruction cross-tab; basis, cutoffTable,
%                   pctlTable, thresholds, recList, reportPctiles
%
%   Bombcell also writes its own cache into <homedir>/bombcell_output.
%
% ---------------------------------------------------------------------
% CURRENT BEHAVIOUR
% ---------------------------------------------------------------------
%   * The percentile gate NEVER REACHES AN ANALYSIS. incl_pctl and keep are
%     descriptive columns; the TSV carries bc_unitType only, and downstream sets
%     pass_bcLabel from that label. The percentile columns describe the shape of
%     the pooled distributions and must not be drawn as a stage of the inclusion
%     funnel.
%   * RPV promotion DOES reach the analysis. A somatic MUA unit is re-classified
%     with the RPV threshold raised to rpvPromoteCutoff and promoted to GOOD iff
%     it flips when ONLY that threshold moves -- i.e. RPV was the sole blocker.
%     NON-SOMA and NOISE cannot flip to somatic GOOD, so they are never promoted.
%     Applied before the TSV is written, so it sets the downstream cohort.
%   * rpvPromoteCutoff = 0.15 is a hand-picked provisional value, not a validated
%     one, and no sensitivity analysis backs it. See the README before relying on
%     the resulting unit counts.
%   * spatialDecaySlope gates 'high': steeper decay = more spatially compact =
%     better. GOOD median slope exceeds NOISE, and corr(slope, SNR) > 0.
%   * Two accountings are kept deliberately separate and must not be conflated --
%     metricPassPctl (percentile, descriptive) and metricPassBC (Bombcell
%     thresholds, gating). In both, nOnlyFail is a metric's MARGINAL cost, the
%     honest number to quote; the plain nFail counts overlap heavily and sum to
%     nothing meaningful.
%   * The label-reconstruction audit rebuilds GOOD from the 9 metricSpec metrics
%     and cross-tabulates against the real label. Bombcell also applies
%     waveform-shape and spikes-missing criteria that are NOT in metricSpec, so
%     exact agreement is not expected; the audit exists to QUANTIFY that residual
%     rather than assume it away.
%   * Panel C and the promotion preview are driven by bcLabelRaw, so they do not
%     move when rpvPromoteCutoff changes. That keeps cutoff selection
%     non-circular: run -> inspect panel C -> set the cutoff -> rerun.
%
% ---------------------------------------------------------------------
% RUN ORDER    THIS -> HKS_cluster_QC_unifiedV7 -> HKS_pool_resultantsV3
% DEPENDS ON   Bombcell (+bc namespace) and its dependencies
% CONVENTIONS  fullfile only, no cd; SETUP -> COMPUTE -> PLOTS; -v7.3
% ---------------------------------------------------------------------
%
% VERSION-DEPENDENT SPOTS, marked in place where they occur:
%   param.splitGoodAndMua_NonSomatic     non-somatic split flag name
%   param.plotGlobal / param.plotDetails figure-suppression field names
%   bc.qm.getQualityUnitType             signature and output shape
%   the RPV classification threshold     field name, see local_relaxRPV
%   the numeric-to-string label map      see local_mapNumericBC
%   the per-metric param field names     see local_bcThresholds

clear; close all

%% =================== BLOCK 0: USER CONFIG ===================

DRIVE = 'W:';                              % drive letter for the data root and the output tree
addpath(genpath(fullfile(DRIVE, 'Kersten', 'Packages_and_Resources')))

% Recordings to include, one row per mouse. There is ONE recording per day,
% so the DATE (YYYY-MM-DD) is the unique key; the "R01" token is NOT a run
% index (it is essentially always R01). Each value is either:
%   'all'                         -> every catgt* recording for the mouse
%   {'2022-04-12','2022-04-13'}   -> only those recording DATES
% Recordings are always returned in chronological (date) order, matching the
% recID = '<date>R0N' convention used by HKS_cluster_QC_unified.
recSpec = { ...
    5, 'all'; ...
    6, 'all'; ...
    7, 'all'; ...
    8, 'all'; ...
    10, 'all' };

ksSubfolder = 'kilosort4';

% Acquisition (Neuropixels 1.0, SpikeGLX defaults)
nChannels       = 385;
ephysSampleRate = 30000;
gainToUV        = 2.3438;           % AP gain 500
region          = 'Cortex';         % unused; the cell-type classifier is steered by
                                    % param.templateDuration_Ctx_threshold below

% Non-somatic split ON so NON-SOMA units get their own label. Field name is
% version-dependent.
splitNonSomatic = true;

% --- Custom RPV promotion ---
% Somatic MUA units that are MUA *solely* because of their RPV fraction are
% promoted to GOOD. Implemented by re-classifying with the RPV classification
% threshold raised to rpvPromoteCutoff (a FRACTION of refractory-period
% violations, same units as fractionRPVs_estimatedTauR) and promoting units that
% flip somatic-MUA -> somatic-GOOD when *only* the RPV threshold changes.
% NON-SOMA and NOISE cannot be promoted this way (their failing metrics are
% untouched), so they stay excluded.
%
% Selecting the cutoff: run, read the GOOD-vs-somatic-MUA split in BLOCK 5
% (panel C) plus the printed promotion preview, set the cutoff, rerun. Panel C
% is driven by the raw labels, so it does not move when the cutoff changes.
% The Bombcell default RPV threshold (printed per recording) is typically ~0.10;
% set rpvPromoteCutoff ABOVE it to promote anything.
enableRPVpromotion = true;
rpvPromoteCutoff   = 0.15;          % provisional, hand-picked; not validated by a
                                    % sweep. Sets the downstream cohort -- see README.

% --- Continuous inclusion gate (descriptive only; does not reach an analysis) ---
inclPctile   = 25;                  % worst 25% per metric fail the gate
reportPctiles = [5 10 25 50 75 90 95];
pctlOnLabelPassOnly = true;         % compute percentiles on non-noise/non-soma
                                    %   units only (recommended). false = pool all.
nanFails     = true;                % NaN on a required metric -> fail that gate

% Per-metric direction: 'high' = larger is better (gate keeps x >= low cut),
%                       'low'  = smaller is better (gate keeps x <= high cut).
% Field names validated against the qMetric table.
metricSpec = { ...
    'nSpikes',                    'high'; ...
    'presenceRatio',              'high'; ...
    'signalToNoiseRatio',         'high'; ...
    'rawAmplitude',               'high'; ...
    'spatialDecaySlope',          'high' ; ...   % steeper = more compact = better
    'fractionRPVs_estimatedTauR', 'low' ; ...
    'nPeaks',                     'low' ; ...
    'nTroughs',                   'low' ; ...
    'maxDriftEstimate',           'low' };

% --- Label-based exclusion (case-insensitive substring match) ---
excludeLabels = {'NOISE', 'NON-SOMA'};   % keep GOOD + (somatic) MUA

forceRerun = false;                 % true = recompute Bombcell even if cached

outDir = fullfile(DRIVE, 'Kersten', 'QC', 'bombcell_pooled');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% =================== BLOCK 1: RESOLVE RECORDING FOLDERS ===================
% Recordings are identified by DATE (one per day). Parse recID ('<date>R0N')
% and recDate ('<date>') from each catgt* folder with the SAME regexes
% HKS_cluster_QC_unified uses, filter by the requested dates (or 'all'), and
% return chronologically. This section only builds the manifest -- inspect
% its printout before running BLOCK 2 (the heavy Bombcell loop).

recList = struct('mouseID', {}, 'recName', {}, 'recID', {}, 'recDate', {}, ...
                 'homedir', {}, 'ksdir', {}, 'rawFile', {}, 'savePath', {});

for i = 1:size(recSpec, 1)
    mouseID = recSpec{i, 1};
    want    = recSpec{i, 2};                 % 'all' | cellstr/string of dates
    ephysDir = fullfile(DRIVE, 'Gulbinaite', 'Resonance_neuropixels', ...
                        'Data_collection', 'Experiment_1', ...
                        sprintf('M%d', mouseID), 'Ephys');

    d = dir(fullfile(ephysDir, 'catgt*'));
    d = d([d.isdir]);
    if isempty(d)
        warning('M%d: no catgt* folders in %s', mouseID, ephysDir); continue;
    end
    recNames = {d.name};
    nD = numel(recNames);

    % --- parse date + recID (project convention) ---
    recDates = repmat("", nD, 1);
    recIDs   = repmat("", nD, 1);
    for k = 1:nD
        tD = regexp(recNames{k}, '(\d{4}-\d{2}-\d{2})',     'tokens', 'once');
        tR = regexp(recNames{k}, '(\d{4}-\d{2}-\d{2}R\d+)', 'tokens', 'once');
        if ~isempty(tD), recDates(k) = string(tD{1}); end
        if ~isempty(tR), recIDs(k) = string(tR{1}); else, recIDs(k) = string(recNames{k}); end
    end
    hasDate = recDates ~= "";
    if any(~hasDate)
        warning('M%d: %d catgt* folder(s) lack a parseable date - ignored.', ...
                mouseID, sum(~hasDate));
    end
    keepIdx = find(hasDate);

    % --- optional date filter ('all' keeps everything) ---
    isAll = (ischar(want) && strcmpi(want, 'all')) || ...
            (isstring(want) && isscalar(want) && strcmpi(want, "all"));
    if ~isAll
        wantDates = string(want);
        sel       = ismember(recDates(keepIdx), wantDates);
        missing   = setdiff(wantDates, recDates(keepIdx));
        if ~isempty(missing)
            warning('M%d: requested dates not found: %s', mouseID, strjoin(missing, ', '));
        end
        keepIdx = keepIdx(sel);
    end

    % --- chronological order (dates are unique per mouse) ---
    [~, ord] = sort(recDates(keepIdx));
    keepIdx  = keepIdx(ord);

    for kk = keepIdx(:)'
        home = fullfile(ephysDir, recNames{kk});
        ap   = dir(fullfile(home, '*ap.bin'));
        if isempty(ap)
            warning('No *ap.bin in %s - skipping.', home); continue;
        end
        e = numel(recList) + 1;
        recList(e).mouseID  = mouseID;
        recList(e).recName  = recNames{kk};
        recList(e).recID    = char(recIDs(kk));
        recList(e).recDate  = char(recDates(kk));
        recList(e).homedir  = home;
        recList(e).ksdir    = fullfile(home, ksSubfolder);
        recList(e).rawFile  = fullfile(home, ap(1).name);
        recList(e).savePath = fullfile(home, 'bombcell_output');
    end
end

nRec = numel(recList);
fprintf('\nResolved %d recordings (chronological):\n', nRec);
for e = 1:nRec
    fprintf('  M%-2d  %-14s  %s\n', recList(e).mouseID, recList(e).recID, recList(e).recName);
end
if nRec == 0, error('No recordings resolved - check DRIVE and recSpec.'); end

%% =================== BLOCK 2: RUN BOMBCELL PER RECORDING ===================
% Collect every unit into one long table keyed by globalID. Also applies the
% custom RPV promotion per recording, at metric-assembly time.

pooled = table();
paramStore = cell(nRec, 1);      % Bombcell's own thresholds, per recording

for e = 1:nRec
    r = recList(e);
    fprintf('\n[%d/%d] M%d  %s  (%s)\n', e, nRec, r.mouseID, r.recID, r.recName);

    if ~exist(r.savePath, 'dir'), mkdir(r.savePath); end

    ephysMetaDir = dir(fullfile(fileparts(r.rawFile), '*.ap.meta'));
    if isempty(ephysMetaDir)
        warning('No .ap.meta in %s - skipping recording.', r.homedir); continue;
    end

    % --- params ---
    param = bc.qm.qualityParamValues(ephysMetaDir, r.rawFile, r.ksdir, [], 4);
    param.gain_to_uV         = gainToUV;
    param.nChannels          = nChannels;
    param.ephys_sample_rate  = ephysSampleRate;
    param.templateDuration_Ctx_threshold = 400;
    param.splitGoodAndMua_NonSomatic = double(splitNonSomatic);  % field name is version-dependent
    param.recomputeQualityMetrics    = double(forceRerun);
    % --- suppress Bombcell's own figures (ours are in BLOCK 5) ---
    for f = {'plotGlobal', 'plotDetails'}                        % field names vary
        if isfield(param, f{1})                                  % across Bombcell versions
            param.(f{1}) = 0;
        elseif e == 1
            warning('param.%s not present in this Bombcell version.', f{1});
        end
    end
    paramStore{e} = param;   % kept for the supplement threshold table

    % --- load KS output ---
    [spikeTimes_samples, spikeClusters, templateWaveforms, templateAmplitudes, ...
     pcFeatures, pcFeatureIdx, channelPositions] = ...
        bc.load.loadEphysData(r.ksdir, r.savePath, 1);

    % --- quality metrics ---
    [qMetric, ~] = bc.qm.runAllQualityMetrics(param, spikeTimes_samples, ...
        spikeClusters, templateWaveforms, templateAmplitudes, ...
        pcFeatures, pcFeatureIdx, channelPositions, r.savePath);

    % alias for classifyCells compatibility
    if ~ismember('waveformDuration_peakTrough_us', qMetric.Properties.VariableNames) && ...
        ismember('waveformDuration_peakTrough',    qMetric.Properties.VariableNames)
        qMetric.waveformDuration_peakTrough_us = qMetric.waveformDuration_peakTrough;
    end

    % --- default string unit-type labels (GOOD / MUA / NOISE / NON-SOMA ...) ---
    bcLabelRaw = local_bcLabels(param, qMetric);    % signature varies; see helper
    Lraw       = upper(bcLabelRaw);

    % --- custom RPV promotion (somatic MUA -> GOOD, RPV-only) ---
    % Re-classify with the RPV threshold raised to rpvPromoteCutoff; a unit is
    % promoted iff it is somatic MUA under the default threshold AND flips to
    % somatic GOOD when *only* the RPV threshold is relaxed (i.e. RPV was the
    % sole blocker). NON-SOMA / NOISE cannot flip to somatic GOOD, so they are
    % never promoted.
    n        = height(qMetric);
    promoted = false(n, 1);
    defRPV   = NaN;
    if enableRPVpromotion
        [paramRPV, rpvFields] = local_relaxRPV(param, rpvPromoteCutoff);
        if isempty(rpvFields)
            warning(['No RPV threshold field found in param for %s. ' ...
                     'RPV promotion SKIPPED (see local_relaxRPV).'], r.recName);
        else
            defRPV       = param.(rpvFields{1});
            bcLabelRelax = local_bcLabels(paramRPV, qMetric);
            Lrelax       = upper(bcLabelRelax);
            isSomMUAraw  = contains(Lraw,   'MUA')  & ~contains(Lraw,   'NON-SOMA');
            isGoodRelax  = contains(Lrelax, 'GOOD') & ~contains(Lrelax, 'NON-SOMA');
            promoted     = isSomMUAraw & isGoodRelax;
        end
    end
    bcLabel           = bcLabelRaw;
    bcLabel(promoted) = "GOOD";

    % --- assemble this recording's rows ---
    phy = qMetric.phy_clusterID(:);
    gid = strings(n, 1);
    for u = 1:n, gid(u) = sprintf('%s__clu%d', r.recName, phy(u)); end

    T = table();
    T.globalID   = gid;
    T.recName    = repmat(string(r.recName), n, 1);
    T.mouseID    = repmat(r.mouseID, n, 1);
    T.clusterID  = qMetric.clusterID(:);
    T.phyID      = phy;
    T.bcLabel    = string(bcLabel(:));       % FINAL (post-promotion) -> TSV + downstream
    T.bcLabelRaw = string(bcLabelRaw(:));    % audit: default Bombcell label
    T.rpvPromoted = promoted(:);             % audit: MUA -> GOOD by relaxed RPV

    % continuous metrics (guard each field so a version mismatch warns, not crashes)
    for m = 1:size(metricSpec, 1)
        f = metricSpec{m, 1};
        if ismember(f, qMetric.Properties.VariableNames)
            T.(f) = qMetric.(f)(:);
        else
            warning('Metric "%s" absent from qMetric (M%d %s) - filling NaN.', ...
                    f, r.mouseID, r.recName);
            T.(f) = nan(n, 1);
        end
    end

    pooled = [pooled; T]; %#ok<AGROW>

    % --- per-recording TSV (FINAL labels) for the downstream unified QC loader ---
    tsv = table(phy, cellstr(bcLabel(:)), ...
                'VariableNames', {'cluster_id', 'bc_unitType'});
    writetable(tsv, fullfile(r.ksdir, 'cluster_bc_unitType.tsv'), ...
               'FileType', 'text', 'Delimiter', '\t');

    fprintf('  %d units | GOOD %d  MUA %d  NOISE %d  NON-SOMA %d', n, ...
        sum(contains(upper(T.bcLabel), 'GOOD') & ~contains(upper(T.bcLabel), 'NON-SOMA')), ...
        sum(contains(upper(T.bcLabel), 'MUA')  & ~contains(upper(T.bcLabel), 'NON-SOMA')), ...
        sum(contains(upper(T.bcLabel), 'NOISE')), ...
        sum(contains(upper(T.bcLabel), 'NON-SOMA')));
    if enableRPVpromotion
        fprintf(' | RPV promote (def %.3g -> %.3g): %d somatic MUA -> GOOD', ...
                defRPV, rpvPromoteCutoff, sum(promoted));
    end
    fprintf('\n');
end

fprintf('\nPooled total: %d units across %d recordings.\n', height(pooled), nRec);
if enableRPVpromotion
    fprintf('RPV promotion (final): %d units promoted somatic MUA -> GOOD (cutoff %.3g).\n', ...
            sum(pooled.rpvPromoted), rpvPromoteCutoff);
end

%% =================== BLOCK 3: LABEL EXCLUSION MASK ===================
% Built on the FINAL (post-promotion) labels, before percentiles so the
% distributions can ignore excluded units.

L = upper(pooled.bcLabel);
excl_label = false(height(pooled), 1);
for j = 1:numel(excludeLabels)
    excl_label = excl_label | contains(L, upper(excludeLabels{j}));
end
pooled.excl_label = excl_label;

fprintf('\nLabel exclusion: %d / %d units flagged (%s)\n', ...
    sum(excl_label), height(pooled), strjoin(excludeLabels, ', '));

%% =================== BLOCK 4: PERCENTILES + GATE ===================
% Distribution/percentiles per metric, then the descriptive inclusion gate.

if pctlOnLabelPassOnly
    basis = ~excl_label;              % percentiles on non-noise/non-soma units
else
    basis = true(height(pooled), 1);  % percentiles on everything
end

nM = size(metricSpec, 1);
metricNames = metricSpec(:, 1);
pctlTable = array2table(nan(nM, numel(reportPctiles)), ...
    'RowNames', metricNames, ...
    'VariableNames', compose('p%d', reportPctiles));
cutoff = nan(nM, 1);                  % the inclusion cutoff per metric

incl_pctl = true(height(pooled), 1);  % passes ALL per-metric gates
passMatPctl = true(height(pooled), nM);   % per-metric, kept for the audit

for m = 1:nM
    f   = metricSpec{m, 1};
    dir = metricSpec{m, 2};
    x   = pooled.(f);
    xb  = x(basis);

    pctlTable{m, :} = prctile(xb, reportPctiles);   % prctile ignores NaN

    switch dir
        case 'high'
            cutoff(m) = prctile(xb, inclPctile);
            pass = x >= cutoff(m);
        case 'low'
            cutoff(m) = prctile(xb, 100 - inclPctile);
            pass = x <= cutoff(m);
        otherwise
            error('metricSpec direction for "%s" must be high|low.', f);
    end

    if nanFails, pass(isnan(x)) = false; else, pass(isnan(x)) = true; end
    passMatPctl(:, m) = pass;
    incl_pctl = incl_pctl & pass;
end

pooled.incl_pctl = incl_pctl;
pooled.keep      = incl_pctl & ~excl_label;

cutoffTable = table(metricNames, metricSpec(:, 2), cutoff, ...
    'VariableNames', {'metric', 'direction', 'cutoff'});

fprintf('\nPercentiles per metric (basis: %s):\n', ...
    ternary(pctlOnLabelPassOnly, 'label-passing units', 'all units'));
disp(pctlTable);
fprintf('Inclusion cutoffs at the %d%% gate:\n', inclPctile);
disp(cutoffTable);
fprintf('Continuous gate passes: %d\n', sum(incl_pctl));
fprintf('KEEP (gate AND not-excluded): %d / %d (%.1f%%)\n', ...
    sum(pooled.keep), height(pooled), 100 * mean(pooled.keep));

%% =================== BLOCK 4b: PER-METRIC ACCOUNTING ===================
% Everything a QC supplement table needs, in two SEPARATE accountings that
% must not be conflated:
%
%   (a) PERCENTILE cutoffs  -> metricPassPctl.  Descriptive only. These numbers
%       describe the shape of the pooled distributions; they exclude nothing
%       downstream (see the CURRENT BEHAVIOUR note in the header).
%   (b) BOMBCELL thresholds -> metricPassBC.    Gating. These are the criteria
%       that actually decide GOOD vs MUA and therefore pass_bcLabel.
%
% "nOnlyFail" = units that fail THAT metric and pass every other metric in the
% same accounting. It is the metric's marginal cost, and it is the honest
% number to quote when a reader asks what a given criterion removed. The plain
% nFail counts are heavily overlapping and will not sum to anything meaningful.

nBasis = sum(basis);

% ---- (a) percentile accounting, restricted to the keep-candidate basis ----
Pb = passMatPctl(basis, :);
nPassP  = sum(Pb, 1)';
nFailP  = nBasis - nPassP;
nOnlyP  = zeros(nM, 1);
for m = 1:nM
    others = true(1, nM); others(m) = false;
    nOnlyP(m) = sum(~Pb(:, m) & all(Pb(:, others), 2));
end
metricPassPctl = table(string(metricNames), string(metricSpec(:, 2)), cutoff, ...
    nPassP, nFailP, 100 * nFailP / max(nBasis, 1), nOnlyP, ...
    'VariableNames', {'metric', 'direction', 'pctlCutoff', ...
                      'nPass', 'nFail', 'pctFail', 'nOnlyFail'});

fprintf('\n--- (a) PERCENTILE accounting (basis = %d keep-candidate units) ---\n', nBasis);
fprintf('    DESCRIPTIVE ONLY: these cutoffs do not propagate downstream.\n');
disp(metricPassPctl);

% ---- (b) Bombcell threshold accounting ----
% Thresholds are read from the stored param. Field names drift across Bombcell
% versions, so every metric is resolved through a CANDIDATE LIST and reported
% as NaN (skipped) if none of the candidates exist. Check the printed table for
% NaNs before quoting it.
[bcThr, bcDir, bcField] = local_bcThresholds(paramStore, metricNames);

passMatBC = true(height(pooled), nM);
for m = 1:nM
    x = pooled.(metricNames{m});
    if isnan(bcThr(m))
        continue                                  % unresolved threshold -> no constraint
    end
    if bcDir(m) == "high", p = x >= bcThr(m); else, p = x <= bcThr(m); end
    p(isnan(x)) = false;                          % same NaN convention as the pctl gate
    passMatBC(:, m) = p;
end

Bb = passMatBC(basis, :);
nPassB = sum(Bb, 1)';
nFailB = nBasis - nPassB;
nOnlyB = zeros(nM, 1);
for m = 1:nM
    others = true(1, nM); others(m) = false;
    nOnlyB(m) = sum(~Bb(:, m) & all(Bb(:, others), 2));
end
metricPassBC = table(string(metricNames), bcDir(:), string(bcField(:)), bcThr(:), ...
    nPassB, nFailB, 100 * nFailB / max(nBasis, 1), nOnlyB, ...
    'VariableNames', {'metric', 'direction', 'paramField', 'bcThreshold', ...
                      'nPass', 'nFail', 'pctFail', 'nOnlyFail'});

fprintf('\n--- (b) BOMBCELL threshold accounting (basis = %d keep-candidate units) ---\n', nBasis);
fprintf('    GATING: these thresholds produce the label that reaches the TSV.\n');
disp(metricPassBC);
if any(isnan(bcThr))
    fprintf(2, '    %d metric(s) had no resolvable param field and were treated as unconstrained:\n      %s\n', ...
        sum(isnan(bcThr)), strjoin(string(metricNames(isnan(bcThr))), ', '));
end

% ---- (c) label reconstruction audit ----
% How much of the GOOD label do the 9 metricSpec metrics account for? Bombcell
% also uses waveform-shape and spikes-missing criteria that are not in
% metricSpec, so a residual is expected. Quantify it rather than assume it.
isGoodReal  = upper(string(pooled.bcLabel)) == "GOOD";
isGoodRecon = all(passMatBC, 2) & ~excl_label;
agree       = mean(isGoodReal(basis) == isGoodRecon(basis));
reconTab    = table( ...
    [sum( isGoodReal(basis) &  isGoodRecon(basis)); sum( isGoodReal(basis) & ~isGoodRecon(basis)); ...
     sum(~isGoodReal(basis) &  isGoodRecon(basis)); sum(~isGoodReal(basis) & ~isGoodRecon(basis))], ...
    'RowNames', {'GOOD & recon-GOOD', 'GOOD & recon-notGOOD', ...
                 'notGOOD & recon-GOOD', 'notGOOD & recon-notGOOD'}, ...
    'VariableNames', {'n'});

fprintf('\n--- (c) Label reconstruction from the 9 metricSpec metrics ---\n');
disp(reconTab);
fprintf('    agreement with the real bcLabel: %.1f%% of the basis (%d units).\n', ...
        100 * agree, nBasis);
fprintf(['    Residual = criteria Bombcell applies that metricSpec omits\n' ...
         '    (waveform shape, %% spikes missing, isolation distance...). Quote the\n' ...
         '    agreement figure in the supplement instead of implying completeness.\n']);

%% =================== BLOCK 5: PLOTS + DISTRIBUTION DIAGNOSTICS ===================
% (A) and (C) are per-metric histograms with the stored inclusion cutoff (red)
% plus 25/50/75 markers (dashed), and the fractionRPVs GOOD-vs-somatic-MUA
% split used to CHOOSE rpvPromoteCutoff.
%
% IMPORTANT: panel (C) and the promotion preview are driven by the RAW labels
% (bcLabelRaw), so they are stable regardless of the current rpvPromoteCutoff.
% That keeps the cutoff-choosing workflow non-circular: run -> inspect (C) ->
% set rpvPromoteCutoff -> rerun.

% ---- (A) per-metric distributions with cutoff + 25/50/75 markers ----
fprintf('\n%-28s %-6s %12s %10s | %10s %10s %10s\n', ...
        'metric','dir','storedCut','%excl@cut','p25','p50','p75');
figure('Color','w','Name','Bombcell pooled metric distributions', ...
       'Position', [80 80 1500 850]);
tl = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
for m = 1:nM
    f   = metricSpec{m, 1};
    dir = metricSpec{m, 2};
    x   = pooled.(f);
    xb  = x(basis & ~isnan(x));
    if isempty(xb), continue; end
    p   = prctile(xb, [25 50 75]);
    if dir=="low", nExcl = sum(xb > cutoff(m)); else, nExcl = sum(xb < cutoff(m)); end
    fprintf('%-28s %-6s %12.4g %9.0f%% | %10.4g %10.4g %10.4g\n', ...
            f, dir, cutoff(m), 100*nExcl/numel(xb), p(1), p(2), p(3));

    nexttile; hold on
    xlo = prctile(xb, 0.5); xhi = prctile(xb, 99.5);
    if xhi <= xlo, xhi = xlo + eps; end
    edges = linspace(xlo, xhi, 41);
    histogram(xb, edges, 'FaceColor', [.6 .6 .6], 'EdgeColor', 'none');
    yl = ylim;
    plot([cutoff(m) cutoff(m)], yl, 'r-', 'LineWidth', 1.5);
    for q = 1:3, plot([p(q) p(q)], yl, 'k--'); end
    xlim([xlo xhi]);
    title(sprintf('%s (%s)\n%.0f%% excl @ cut', f, dir, 100*nExcl/numel(xb)), ...
          'Interpreter', 'none', 'FontSize', 8);
    xlabel(f, 'Interpreter', 'none'); ylabel('# units'); box off
end
title(tl, sprintf('Pooled distributions (basis: %s), %d%% cutoff red, 25/50/75 dashed', ...
    ternary(pctlOnLabelPassOnly,'keep candidates','all'), inclPctile));

% ---- (B) unit-type composition per recording (FINAL labels) ----
figure('Color','w','Name','Bombcell labels per recording', 'Position', [120 120 900 450]);
recNames = unique(pooled.recName, 'stable');
cats = {'GOOD', 'MUA', 'NOISE', 'NON-SOMA'};
counts = zeros(numel(recNames), numel(cats));
for a = 1:numel(recNames)
    Lr = upper(pooled.bcLabel(pooled.recName == recNames(a)));
    counts(a, 1) = sum(contains(Lr, 'GOOD') & ~contains(Lr, 'NON-SOMA'));
    counts(a, 2) = sum(contains(Lr, 'MUA')  & ~contains(Lr, 'NON-SOMA'));
    counts(a, 3) = sum(contains(Lr, 'NOISE'));
    counts(a, 4) = sum(contains(Lr, 'NON-SOMA'));
end
bar(counts, 'stacked');
set(gca, 'XTick', 1:numel(recNames), 'XTickLabel', recNames, ...
    'TickLabelInterpreter', 'none'); xtickangle(30);
ylabel('# units'); legend(cats, 'Location', 'eastoutside');
title('Unit-type composition per recording (post-promotion)');

% ---- (C) fractionRPVs: GOOD vs somatic-MUA (choose rpvPromoteCutoff) ----
% RAW labels, so this figure/preview does not move when rpvPromoteCutoff changes.
rpvName = 'fractionRPVs_estimatedTauR';
if ismember(rpvName, pooled.Properties.VariableNames)
    Lr0       = upper(pooled.bcLabelRaw);
    isNonSoma = contains(Lr0, 'NON-SOMA');
    isNoise   = contains(Lr0, 'NOISE');
    isGoodRaw = contains(Lr0, 'GOOD') & ~isNonSoma;
    isMUAsom  = contains(Lr0, 'MUA')  & ~isNonSoma;      % somatic MUA only
    rpv       = double(pooled.(rpvName));
    basisR    = ~isNoise & ~isNonSoma;

    figure('Color','w','Name','fractionRPVs: GOOD vs somatic-MUA (RPV promotion cutoff)');
    hold on
    xhi = prctile(rpv(basisR & ~isnan(rpv)), 99);
    edges = linspace(0, max(xhi, eps), 51);
    histogram(rpv(isGoodRaw & ~isnan(rpv)), edges, 'FaceColor', [.15 .5 .15], 'EdgeColor','none', 'FaceAlpha',.6);
    histogram(rpv(isMUAsom  & ~isnan(rpv)), edges, 'FaceColor', [.75 .35 .1], 'EdgeColor','none', 'FaceAlpha',.6);
    xline(0.05, 'k--'); xline(0.10, 'k:');               % common RPV cutoffs, for orientation
    xline(rpvPromoteCutoff, 'r-', 'LineWidth', 1.5);     % the ACTIVE promotion cutoff
    xlabel('fractionRPVs (estimated tauR)'); ylabel('# units'); box off
    legend({'GOOD (raw)','somatic MUA (raw)','0.05','0.10', ...
            sprintf('active cutoff %.3g', rpvPromoteCutoff)}, 'Location','northeast');
    title('somatic MUA below the cutoff (RPV-only) -> promoted to GOOD');

    fprintf('\nMUA promotion preview (raw somatic MUA with fractionRPVs <= cutoff):\n');
    for c = [0.02 0.05 0.08 0.10 0.15]
        fprintf('  cutoff %.2f -> %d of %d somatic-MUA units below RPV cut\n', ...
                c, sum(isMUAsom & rpv <= c), sum(isMUAsom));
    end
    fprintf(['NOTE: this preview counts raw somatic MUA below the RPV cut; the actual\n' ...
             '      promotion also requires the unit to pass every OTHER GOOD criterion\n' ...
             '      (enforced by re-classification), so realised promotions <= preview.\n']);
end

%% =================== BLOCK 6: SAVE ===================

thresholds = struct('inclPctile', inclPctile, 'excludeLabels', {excludeLabels}, ...
    'metricSpec', {metricSpec}, 'pctlOnLabelPassOnly', pctlOnLabelPassOnly, ...
    'nanFails', nanFails, 'splitNonSomatic', splitNonSomatic, ...
    'enableRPVpromotion', enableRPVpromotion, 'rpvPromoteCutoff', rpvPromoteCutoff);

outFile = fullfile(outDir, 'bombcell_pooled_QC.mat');
save(outFile, 'pooled', 'pctlTable', 'cutoffTable', 'thresholds', ...
     'recList', 'reportPctiles', ...
     'paramStore', 'metricPassPctl', 'metricPassBC', ...
     'passMatPctl', 'passMatBC', 'reconTab', 'basis', '-v7.3');

fprintf('\nSaved:\n  %s\n', outFile);
fprintf('  (+ per-recording cluster_bc_unitType.tsv with FINAL labels in each kilosort4 folder)\n');

%% =================== LOCAL FUNCTIONS ===================

function [thr, dirn, fieldUsed] = local_bcThresholds(paramStore, metricNames)
% Resolve Bombcell's OWN classification threshold for each metric in
% metricNames, out of the stored param structs.
%
% param field names drift across Bombcell versions, so each metric is resolved
% through an ordered CANDIDATE LIST and left NaN if none is present. Anything
% NaN is reported as unconstrained rather than silently defaulted.
%
% Also checks that the thresholds are identical across recordings and WARNS if
% they are not, using the first value. They should agree (param is rebuilt from
% the same defaults each time); if a recording ever picks up different values
% the pooled per-metric counts are not interpretable.

    % metric -> {direction, candidate param fields (first hit wins)}
    map = { ...
        'nSpikes',                    'high', {'minNumSpikes','min_nSpikes'}; ...
        'presenceRatio',              'high', {'minPresenceRatio','min_presenceRatio'}; ...
        'signalToNoiseRatio',         'high', {'minSNR','minSignalToNoiseRatio','min_SNR'}; ...
        'rawAmplitude',               'high', {'minAmplitude','min_rawAmplitude'}; ...
        'spatialDecaySlope',          'high', {'minSpatialDecaySlope','minSpatialDecaySlopeExp', ...
                                               'min_spatialDecaySlope'}; ...
        'fractionRPVs_estimatedTauR', 'low',  {'maxRPVviolations','maxRPV','max_fractionRPVs'}; ...
        'nPeaks',                     'low',  {'maxNPeaks','max_nPeaks'}; ...
        'nTroughs',                   'low',  {'maxNTroughs','max_nTroughs'}; ...
        'maxDriftEstimate',           'low',  {'maxDrift','maxDriftEstimate','max_drift'} };

    nM        = numel(metricNames);
    thr       = nan(nM, 1);
    dirn      = strings(nM, 1);
    fieldUsed = repmat("<unresolved>", nM, 1);

    P = paramStore(~cellfun(@isempty, paramStore));
    if isempty(P)
        warning('local_bcThresholds: paramStore is empty (no recording ran) - all NaN.');
        dirn(:) = "high"; return
    end

    for m = 1:nM
        row = find(strcmp(map(:, 1), metricNames{m}), 1);
        if isempty(row)
            warning('local_bcThresholds: no param mapping for "%s".', metricNames{m});
            dirn(m) = "high"; continue
        end
        dirn(m) = string(map{row, 2});
        cands   = map{row, 3};

        vals = nan(numel(P), 1);
        hit  = '';
        for c = 1:numel(cands)
            if isfield(P{1}, cands{c})
                hit = cands{c};
                for k = 1:numel(P)
                    if isfield(P{k}, hit) && isscalar(P{k}.(hit)) && isnumeric(P{k}.(hit))
                        vals(k) = double(P{k}.(hit));
                    end
                end
                break
            end
        end
        if isempty(hit), continue; end

        fieldUsed(m) = string(hit);
        u = unique(vals(~isnan(vals)));
        if numel(u) > 1
            warning(['local_bcThresholds: param.%s differs across recordings ' ...
                     '(%s). Using the first; the pooled per-metric counts are ' ...
                     'only interpretable if these agree.'], hit, ...
                     strjoin(compose('%g', u(:)'), ', '));
        end
        if ~isempty(u), thr(m) = u(1); end
    end
end

function s = local_bcLabels(param, qMetric)
% Return per-unit STRING labels (GOOD / MUA / NOISE / NON-SOMA ...).
% getQualityUnitType has changed shape across Bombcell versions. Primary path
% uses the string output; falls back to a numeric->string map.
    n = height(qMetric);
    s = strings(n, 1);
    got = false;

    % newer: [numericType, stringType] = bc.qm.getQualityUnitType(param, qMetric)
    try
        [~, str] = bc.qm.getQualityUnitType(param, qMetric);
        s = string(str(:)); got = true;
    catch
    end

    % single-output form (string or numeric)
    if ~got
        try
            out = bc.qm.getQualityUnitType(param, qMetric);
            if isstring(out) || iscellstr(out) || ischar(out)
                s = string(out(:));
            else
                s = local_mapNumericBC(double(out(:)));   % check the code map against your build
            end
            got = true;
        catch
        end
    end

    if ~got
        error(['Could not obtain Bombcell string labels. Check the ' ...
               'bc.qm.getQualityUnitType signature for your version.']);
    end
    if numel(s) ~= n
        error('Label count (%d) does not match unit count (%d).', numel(s), n);
    end
end

function [paramOut, hitFields] = local_relaxRPV(param, cutoff)
% Raise the Bombcell RPV *classification* threshold to `cutoff` (a fraction of
% refractory-period violations, same scale as fractionRPVs_estimatedTauR) so
% getQualityUnitType re-thresholds only RPV.
% The threshold field name has drifted across Bombcell builds; set every
% candidate field that actually exists in `param`. Returns the names hit so the
% caller can (a) read the pre-override default and (b) warn on a no-op.
% Canonical field is 'maxRPVviolations'; adjust the list to your build if the
% per-recording print shows defRPV = NaN.
    candidates = {'maxRPVviolations', 'maxRPVViolations', ...
                  'maxfractionRPVs', 'maxFractionRPVs', 'RPVthreshold'};
    paramOut  = param;
    hitFields = {};
    for i = 1:numel(candidates)
        f = candidates{i};
        if isfield(paramOut, f)
            paramOut.(f) = cutoff;
            hitFields{end+1} = f; %#ok<AGROW>
        end
    end
end

function s = local_mapNumericBC(code)
% Numeric -> string fallback. Check against your Bombcell: with the
% non-somatic split ON the common coding is
%   0 NOISE | 1 GOOD | 2 MUA | 3 NON-SOMA (good) | 4 NON-SOMA MUA
    map = ["NOISE"; "GOOD"; "MUA"; "NON-SOMA"; "NON-SOMA MUA"];
    s = strings(numel(code), 1);
    for i = 1:numel(code)
        c = code(i);
        if c >= 0 && c <= numel(map) - 1, s(i) = map(c + 1);
        else, s(i) = "UNKNOWN"; end
    end
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
