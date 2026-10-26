%% HKS_unit_rastersV3.m
% Frequency-ordered flicker rasters for a hand-picked list of units (specified
% by globalID). One unit per region (V1 / dLGN / CA1), drawn as a combined
% 1 x nUnit composite panel.
%
% CA1 is a true null -- no significant following at any frequency -- so its
% raster reads as ABSENCE of frequency structure. That presence-versus-absence
% contrast is the reason to include it.
%
% Same raster logic as HKS_zeta_based_rastersV6 (trials sorted by stimulus
% frequency, y-axis relabelled by frequency instead of trial index), but:
%   - units are named explicitly instead of "top-5 by ZETA";
%   - spikes and stim times come from <recName>_cluster_QC.mat (AP clock,
%     already T0- and rate-corrected), so no NI/CatGT re-sync is needed;
%   - alternating frequency bands + optional cycle ticks for readability.
%
% ---------------------------------------------------------------------
% INPUTS
% ---------------------------------------------------------------------
%   cfg.targets    globalID strings, "<recName>__clu<id>". recName and mouseID
%                  are parsed back out of the string, so the folder must exist
%                  at <dataRoot>/M<n>/Ephys/<recName>.
%   <homedir>/<recName>_cluster_QC.mat
%                  clusterTable, spikes (AP clock) and stim. Falls back to a
%                  bare cluster_QC.mat. sQC.stim is required -- this script does
%                  not reparse structEP.
%
% ---------------------------------------------------------------------
% OUTPUT  ->  <DRIVE>\Kersten\Figures\unit_rasters\
% ---------------------------------------------------------------------
%   raster_panel_<win>.png / .fig / .pdf          the composite (default)
%   raster_<globalID>_<win>.png / .fig / .pdf     per unit, if cfg.doIndividual
%
% ---------------------------------------------------------------------
% CURRENT BEHAVIOUR
% ---------------------------------------------------------------------
%   * The composite is the deliverable; per-unit figures and the short 500 ms
%     window are behind toggles (cfg.doIndividual, cfg.windows) rather than
%     deleted, so either comes back without editing [P1].
%   * An unresolvable target is SKIPPED with a warning rather than aborting the
%     run, so the remaining panels still render while a target is being settled.
%     nUnit is recounted from the resolved set before plotting.
%   * The frequency axis is capped by cfg.dblFreqHigh, well below the 243 Hz
%     acquisition band: above roughly 50 Hz the bands are too thin to read as a
%     raster. The full band is still analysed everywhere else.
%   * Panel titles are one line -- region (+ V1 layer), mouse and cluster.
%     Depth-from-surface and the BC label are QC provenance rather than figure
%     content, so they are left out; spike count is behind cfg.titleShowSpikes.
%   * Type is sized for print, with cfg.fontSize, cfg.labelFontSize and
%     cfg.titleFontSize tunable independently.
%   * Trials are sorted by frequency ascending, chronological within frequency,
%     and the y-axis is relabelled with the frequency of each band rather than
%     trial index. Row 1 is the lowest frequency unless cfg.lowFreqAtTop is off.
%   * Cycle ticks are drawn only up to cfg.cycleTickMaxFreq; above that the
%     lines merge and obscure the spikes they are meant to contextualise.
%   * An empty second element in a cfg.windows entry means "run to stimulus
%     offset", taken per unit from the median trial duration.
%
% ---------------------------------------------------------------------
% RUN ORDER    HKS_cluster_QC_unifiedV7 -> THIS   (reads QC files directly, not
%              the pooled resultant file)
% DEPENDS ON   nothing beyond base MATLAB
% ---------------------------------------------------------------------

clear; close all;

%% ####################################################################
%% ##  [S1] SETUP
%% ####################################################################
cfg = struct();

% --- units to plot: globalID = "<recName>__clu<id>" ------------------
% One representative unit per region, chosen by HKS_pick_exemplarsV1 rather than
% by hand (see README). dLGN and V1 are top-ranked followers; CA1 is picked as
% TYPICAL, not top-ranked, since there is no following in CA1 to rank on. V1 and
% CA1 come from the same recording by design, so the presence-versus-absence
% contrast is not confounded by session.
 cfg.targets = { ...
    'catgt_RecM7_2022-04-14R01_g0__clu473', ...   % V1
    'catgt_RecM8_2022-04-14R01_g0__clu216', ...   % dLGN
    'catgt_RecM7_2022-04-14R01_g0__clu353'  ...   % CA1
    };
% --- paths ----------------------------------------------------------
cfg.DRIVE       = 'W:';            % drive letter for the data root and the output tree
cfg.dataRoot    = fullfile(cfg.DRIVE, 'Gulbinaite', 'Resonance_neuropixels', ...
                           'Data_collection', 'Experiment_1');
cfg.figDir      = fullfile(cfg.DRIVE, 'Kersten', 'Figures', 'unit_rasters');
cfg.saveFigs    = true;

% --- stimulus / trial selection -------------------------------------
cfg.dblFreqLow  = 2;               % Hz, flicker band kept
cfg.dblFreqHigh = 40;              % Hz, upper cap for legibility (see header)
% Curated subset for a clean, evenly labelled y-axis. Each value must be a
% frequency that was actually presented, otherwise it is simply absent --
% [C1] warns per unit listing any requested-but-missing frequency, so if
% 50 Hz was not a stimulus, swap it for the nearest presented value.
cfg.freqSubset  = [];   % [] = all trials in [low,high]

% --- time windows (one figure per window) ----------------------------
% [] as the 2nd element of a window = "run to stimulus offset" (full trial).
% Add { [-0.05 0.5], ... } / 'first500ms' here to get the short window back;
% no change is needed in [P1].
cfg.windows     = { [-0.1 2] };
cfg.winNames    = { 'firsthalf' };

% --- what to render --------------------------------------------------
cfg.doIndividual = false;   % true: also write one figure per unit
cfg.doComposite  = true;    % the 1 x nUnit presentation panel

% --- raster appearance ----------------------------------------------
cfg.lowFreqAtTop    = true;    % true: lowest freq on top, highest at bottom
cfg.tickHeight      = 0.85;    % spike tick height, in trial units
cfg.tickWidth       = 0.5;
cfg.tickColor       = [0 0 0];
cfg.bandShade       = [0.93 0.93 0.96];  % alternating frequency-band shading
cfg.showCycleTicks  = true;    % faint vertical lines at stimulus cycle onsets
cfg.cycleTickMaxFreq= 20;      % ...only for frequencies <= this (else unreadable)
cfg.cycleTickColor  = [0.85 0.35 0.35];
cfg.labelEveryFreq  = true;    % false: label every 2nd frequency band
                               % NOTE: with print-size type and ~20 bands
                               % in 2-40 Hz the y labels can collide. If they
                               % do, set this false before shrinking fontSize.
cfg.fontSize        = 13;      % tick labels
cfg.labelFontSize   = 15;      % x/y axis labels
cfg.titleFontSize   = 16;      % panel titles
cfg.titleShowSpikes = false;   % append "| N spikes" to the panel title

if cfg.saveFigs && ~exist(cfg.figDir, 'dir'), mkdir(cfg.figDir); end

%% ####################################################################
%% ##  [C1] COMPUTE : gather spikes + trials for each target unit
%% ####################################################################
nUnit = numel(cfg.targets);
U     = struct([]);                       % per-unit bundle for the plot stage
kUnit = 0;                                % # of successfully resolved units

for u = 1:nUnit
    gid = cfg.targets{u};
    try
        % --- parse globalID -> recName + clusterID ------------------
        tok = regexp(gid, '^(?<rec>.+)__clu(?<clu>\d+)$', 'names');
        if isempty(tok), error('Cannot parse globalID "%s".', gid); end
        recName   = tok.rec;
        clusterID = str2double(tok.clu);

        % --- mouse from recName (…RecM7_… -> 7) ---------------------
        mtok = regexp(recName, 'RecM(?<m>\d+)', 'names');
        if isempty(mtok), error('Cannot infer mouseID from recName "%s".', recName); end
        mouseID  = str2double(mtok.m);
        ephysDir = fullfile(cfg.dataRoot, sprintf('M%d', mouseID), 'Ephys');

        % --- resolve recording folder -------------------------------
        homedir = fullfile(ephysDir, recName);
        if ~exist(homedir, 'dir')
            error('No folder %s under %s.', recName, ephysDir);
        end

        % --- QC bundle ----------------------------------------------
        qcFile = fullfile(homedir, sprintf('%s_cluster_QC.mat', recName));
        if ~exist(qcFile, 'file')
            qcFile = fullfile(homedir, 'cluster_QC.mat');
        end
        if ~exist(qcFile, 'file')
            error('No *_cluster_QC.mat in %s.', homedir);
        end
        sQC = load(qcFile);
        C   = sQC.clusterTable;

        row = find(C.clusterID == clusterID, 1);
        if isempty(row)
            error('Cluster %d not present in %s.', clusterID, qcFile);
        end

        % --- spikes (AP clock; already rate- and T0-corrected) ------
        S   = sQC.spikes;
        spk = S.spikeTime(S.clusterID == clusterID);
        spk = sort(spk(:));

        % --- flicker trials (AP clock) ------------------------------
        if ~(isfield(sQC, 'stim') && isfield(sQC.stim, 'flickerOn_AP'))
            error(['sQC.stim missing in %s. Re-run HKS_cluster_QC_unifiedV7 ' ...
                   '(the sQC.stim = stim save).'], qcFile);
        end
        onT  = sQC.stim.flickerOn_AP(:);
        frex = sQC.stim.flickerFreq(:);
        if isfield(sQC.stim, 'flickerOff_AP')
            offT = sQC.stim.flickerOff_AP(:);
        else
            offT = onT + 4.988;
        end

        keep = frex >= cfg.dblFreqLow & frex <= cfg.dblFreqHigh & ~isnan(frex);
        if ~isempty(cfg.freqSubset)
            keep = keep & ismember(frex, cfg.freqSubset(:));
        end
        onT = onT(keep);  offT = offT(keep);  frex = frex(keep);

        % --- warn if any requested subset frequency was not presented
        if ~isempty(cfg.freqSubset)
            want    = cfg.freqSubset(cfg.freqSubset >= cfg.dblFreqLow & ...
                                     cfg.freqSubset <= cfg.dblFreqHigh);
            missing = setdiff(want(:), unique(frex));
            if ~isempty(missing)
                warning('%s: requested freq(s) not presented, dropped: %s Hz', ...
                        gid, mat2str(missing(:)'));
            end
        end
        if isempty(onT)
            error('No trials left for %s after band/subset filtering.', gid);
        end

        % --- trial order: frequency ascending, chronological within -
        [~, ord] = sortrows([frex, onT]);
        onT = onT(ord);  offT = offT(ord);  frex = frex(ord);

        [freqList, ~, grpIdx] = unique(frex);  % grpIdx: trial -> freq band

        % --- metadata for the title ---------------------------------
        meta = struct();
        meta.globalID  = gid;
        meta.recName   = recName;
        meta.mouseID   = mouseID;
        meta.clusterID = clusterID;
        meta.region    = getCol(C, row, 'region',       "?");
        meta.v1Layer   = getCol(C, row, 'v1Layer',      "");
        meta.depth     = getCol(C, row, 'depthFromSurface', NaN);
        meta.bcType    = getCol(C, row, 'bcType',       "");
        meta.zetaP     = getCol(C, row, 'zetaP_pulse',  NaN);
        meta.nSpikes   = numel(spk);
        meta.stimDur   = median(offT - onT);

        kUnit = kUnit + 1;
        U(kUnit).spk      = spk;   %#ok<*SAGROW>
        U(kUnit).onT      = onT;
        U(kUnit).offT     = offT;
        U(kUnit).frex     = frex;
        U(kUnit).freqList = freqList;
        U(kUnit).grpIdx   = grpIdx;
        U(kUnit).meta     = meta;

        fprintf('%-42s  %-4s  %5d spikes  %3d trials  %2d freqs  (%.3f s stim)\n', ...
            gid, meta.region, meta.nSpikes, numel(onT), numel(freqList), meta.stimDur);
    catch ME
        warning('Skipping target "%s": %s', gid, ME.message);
    end
end

nUnit = numel(U);              % from here on, only successfully resolved units
if nUnit == 0
    error('No target units could be resolved -- check cfg.targets / paths.');
end

%% ####################################################################
%% ##  [P1] PLOTS : one figure per unit, plus a combined panel
%% ####################################################################
for w = 1:numel(cfg.windows)
    win     = cfg.windows{w};
    winName = cfg.winNames{w};

    % --- individual figures (off by default in V3) --------------------
    if cfg.doIndividual
        for u = 1:nUnit
            figure('Color', 'w', 'Units', 'normalized', ...
                   'Position', [0.08 + 0.06*u, 0.12, 0.34, 0.72]);
            ax = axes(); %#ok<LAXES>
            drawFreqRaster(ax, U(u), win, cfg);

            if cfg.saveFigs
                base = sprintf('raster_%s_%s', U(u).meta.globalID, winName);
                savefig(gcf, fullfile(cfg.figDir, [base '.fig']));
                exportgraphics(gcf, fullfile(cfg.figDir, [base '.png']), 'Resolution', 300);
                exportgraphics(gcf, fullfile(cfg.figDir, [base '.pdf']), 'ContentType', 'vector');
            end
        end
    end

    % --- combined presentation panel ---------------------------------
    if cfg.doComposite
        % taller and wider than V2: the larger type needs the room, and the
        % y axis carries ~20 frequency labels.
        panelW = min(0.96, 0.26*nUnit + 0.14);   % widen with panel count
        figure('Color', 'w', 'Units', 'normalized', 'Position', [0.03 0.06 panelW 0.86]);
        % Only blank the repeated y tick labels if every unit really was shown
        % the same frequency set. Units come from different recordings, so a
        % mismatched set must keep its own labels or the panel lies.
        sameFreq = all(arrayfun(@(k) isequal(U(k).freqList, U(1).freqList), 1:nUnit));
        if ~sameFreq
            warning(['Units differ in presented frequencies -- per-panel y ' ...
                     'labels kept. Check the band before using as one figure.']);
        end

        tl = tiledlayout(1, nUnit, 'TileSpacing', 'compact', 'Padding', 'compact');
        for u = 1:nUnit
            ax = nexttile(tl);
            drawFreqRaster(ax, U(u), win, cfg);
            if u > 1
                ylabel(ax, '');                       % shared frequency axis
                if sameFreq, set(ax, 'YTickLabel', []); end
            end
        end

        if cfg.saveFigs
            base = sprintf('raster_panel_%s', winName);
            savefig(gcf, fullfile(cfg.figDir, [base '.fig']));
            exportgraphics(gcf, fullfile(cfg.figDir, [base '.png']), 'Resolution', 300);
            exportgraphics(gcf, fullfile(cfg.figDir, [base '.pdf']), 'ContentType', 'vector');
        end
    end
end

fprintf('\nDone. Figures in %s\n', cfg.figDir);

%% ####################################################################
%% ##  LOCAL FUNCTIONS
%% ####################################################################
function drawFreqRaster(ax, Uu, win, cfg)
% Frequency-ordered raster on axes AX for unit bundle Uu, window WIN.
% Trials are already sorted low->high frequency; row 1 = lowest frequency.

spk      = Uu.spk;
onT      = Uu.onT;
frex     = Uu.frex;
freqList = Uu.freqList;
grpIdx   = Uu.grpIdx;
M        = Uu.meta;

nTrial = numel(onT);
nFreq  = numel(freqList);

% window: [] as post = run to stimulus offset
tPre  = win(1);
if numel(win) < 2 || isempty(win(2)), tPost = M.stimDur;
else,                                 tPost = win(2); end

hold(ax, 'on');

% --- alternating frequency-band shading ------------------------------
for f = 1:2:nFreq
    rows = find(grpIdx == f);
    patch(ax, [tPre tPost tPost tPre], ...
              [min(rows)-0.5, min(rows)-0.5, max(rows)+0.5, max(rows)+0.5], ...
              cfg.bandShade, 'EdgeColor', 'none');
end

% --- stimulus cycle onsets (low frequencies only) ---------------------
if cfg.showCycleTicks
    for f = 1:nFreq
        fHz = freqList(f);
        if fHz > cfg.cycleTickMaxFreq, continue; end
        rows  = find(grpIdx == f);
        y0    = min(rows) - 0.5;  y1 = max(rows) + 0.5;
        cyc   = 0 : 1/fHz : tPost;
        cyc   = cyc(cyc >= tPre & cyc <= tPost);
        Xc = [repmat(cyc, 2, 1); nan(1, numel(cyc))];
        Yc = [repmat([y0; y1], 1, numel(cyc)); nan(1, numel(cyc))];
        plot(ax, Xc(:), Yc(:), '-', 'Color', [cfg.cycleTickColor 0.35], 'LineWidth', 0.4);
    end
end

% --- spike ticks -------------------------------------------------------
% For each trial, spikes falling in [onset+tPre, onset+tPost], expressed
% relative to onset. Built as NaN-separated segments -> one plot call.
X = []; Y = [];
h = cfg.tickHeight / 2;
for t = 1:nTrial
    s = spk(spk >= onT(t) + tPre & spk <= onT(t) + tPost) - onT(t);
    if isempty(s), continue; end
    s = s(:)';
    X = [X, [s; s; nan(1, numel(s))]];              %#ok<AGROW>
    Y = [Y, [repmat([t-h; t+h], 1, numel(s)); nan(1, numel(s))]]; %#ok<AGROW>
end
if ~isempty(X)
    plot(ax, X(:), Y(:), '-', 'Color', cfg.tickColor, 'LineWidth', cfg.tickWidth);
end

% --- stimulus onset (and offset, if in view) --------------------------
xline(ax, 0, 'k-', 'LineWidth', 0.8);
if tPost >= M.stimDur - 1e-6
    xline(ax, M.stimDur, 'k:', 'LineWidth', 0.8);
end

% --- axes: y ticks at the CENTRE of each frequency band ---------------
ctr = arrayfun(@(f) mean(find(grpIdx == f)), 1:nFreq);
lbl = arrayfun(@(f) num2str(freqList(f)), 1:nFreq, 'UniformOutput', false);
if ~cfg.labelEveryFreq
    lbl(2:2:end) = {''};
end
set(ax, 'YTick', ctr, 'YTickLabel', lbl, 'FontSize', cfg.fontSize, ...
        'TickDir', 'out', 'Box', 'off', 'Layer', 'top');
ylim(ax, [0.5, nTrial + 0.5]);
xlim(ax, [tPre, tPost]);
if cfg.lowFreqAtTop, set(ax, 'YDir', 'reverse'); end

xlabel(ax, 'Time from flicker onset (s)', 'FontSize', cfg.labelFontSize);
ylabel(ax, 'Stimulus frequency (Hz)',     'FontSize', cfg.labelFontSize);

% --- title -------------------------------------------------------------
% One line. Depth-from-surface and the BC/KS label are QC provenance and
% live in the methods / cluster table, not on the figure. Spike count is
% optional (cfg.titleShowSpikes).
regStr = char(M.region);
if strcmp(regStr, 'V1') && strlength(string(M.v1Layer)) > 0
    regStr = sprintf('V1 %s', char(M.v1Layer));
end
tStr = sprintf('%s  |  M%d clu%d', regStr, M.mouseID, M.clusterID);
if cfg.titleShowSpikes
    tStr = sprintf('%s  |  %d spikes', tStr, M.nSpikes);
end
title(ax, tStr, 'FontSize', cfg.titleFontSize, 'FontWeight', 'normal');
end

function v = getCol(T, row, name, dflt)
% Fetch T.(name)(row) if the column exists, else DFLT.
if ismember(name, T.Properties.VariableNames)
    v = T.(name)(row);
    if iscell(v), v = v{1}; end
else
    v = dflt;
end
end
