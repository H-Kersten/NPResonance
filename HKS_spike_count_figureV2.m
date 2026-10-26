%% HKS_spike_count_figureV2.m
% Raw spike-count / exposure figure from the pooled resultant file.
% Downstream consumer of HKS_pool_resultantsV3.m -- NO SPIKES TOUCHED.
% Reads ONLY Nc (spike count per cluster x frequency), Nn (the blank-side count)
% and Tc (observation duration per cluster x frequency). No Zc, no metrics, no
% zetatest.
%
% ---------------------------------------------------------------------
% WHY THIS FIGURE EXISTS (three jobs, one panel each)
% ---------------------------------------------------------------------
%  (1) YIELD / SANITY. How many spikes does each cluster actually contribute,
%      per area? Notably, CA1 sits in the same yield range as dLGN/V1, so the
%      CA1 null cannot be dismissed as a power problem. [P-1 A, P-1 B]
%  (2) THE VARIABLE DENOMINATOR. The minSpkPerBin guard drops a DIFFERENT set
%      of cells at every frequency, so the "n" behind every pooled PPC point
%      differs. This is the honest disclosure of why a SEM-across-cells band is
%      anticonservative as an animal-level statement (the animal-level object
%      remains the hierarchical bootstrap in the descriptor script).
%      It is ALSO the counting-floor defence of f_cutoff: the guard can only
%      sculpt a cutoff if it thins cells preferentially NEAR that frequency, so
%      a denominator that stays flat with frequency means an observed cutoff is
%      fidelity loss, not attrition. [P-1 C]
%  (3) THE |F1| RATE CONFOUND, MADE VISIBLE. spikes/cycle = Nc ./ (Tc .* f) =
%      rate/f has a built-in ~1/f envelope. |F1| in per-cycle units is bounded
%      by this quantity, so |F1| inherits a 1/f roll-off that has nothing to do
%      with genuine following. This is the quantitative argument for preferring
%      PPC (rate-independent, null-centred) over |F1| for any frequency-
%      dependence claim. It does NOT constrain PPC, which pools phase over all
%      cycles and is untroubled by cells that fire on a minority of them
%      (see the cycle-skipping / f/2 discussion). [P-1 D]
%
% ---------------------------------------------------------------------
% CURRENT BEHAVIOUR
% ---------------------------------------------------------------------
%   * TWO KINDS OF MISSING, drawn differently in [P-1 A]:
%       NaN Nc            -> that frequency was not presented in that recording
%                            (union frequency grid) -> WHITE (AlphaData = 0)
%       Nc < minSpkPerBin -> presented, but below the guard -> GREY OVERLAY
%     Collapsing these two would misrepresent both yield and exclusion.
%   * The guard shown is the PAIRED guard (Nc AND Nn >= minSpkPerBin), which is
%     the denominator the descriptor f_cutoff test actually uses. Nn is the
%     binding side and drops roughly 40% of bins at high frequency, so the
%     paired guard, not the real-only guard, is what the survival panel must
%     show. Set showPairedGuard = false only to inspect the real-trial guard in
%     isolation.
%   * minSpkPerBin here must match the value used by the descriptor and guard-
%     sweep scripts, or panels A and C document a guard that was never applied.
%     Single knob in [S4].
%   * Panel A uses a COLUMN-INDEX x-axis (1..nFreq), not a true log-Hz axis,
%     because imagesc cannot render non-uniform column spacing. The frequency
%     grid is near-log-spaced, so the index axis is visually close to log-Hz;
%     ticks are labelled with the true frequencies. Panels C and D use a genuine
%     log-Hz axis with shared ticks. Do not compare slopes or x-positions across
%     the two axis types (A vs C/D).
%   * The inclusion mask can only FURTHER restrict what the pool kept, since the
%     pool is already QC-pre-filtered. Same ceiling caveat as the descriptor
%     script.
%   * Panel D medians are taken over ALL presented cells, NOT guard-gated: the
%     point of the panel is to show the exposure the guard is reacting to, so
%     gating there would be circular.
%
% ---------------------------------------------------------------------
% INPUT   <DRIVE>\Kersten\Pooled\resonance_resultants_pool.mat
% OUTPUT  <DRIVE>\Kersten\Figures\spike_counts_exposure.png / .fig
% RUN ORDER    HKS_pool_resultantsV3 -> THIS
% DEPENDS ON   Statistics Toolbox (prctile). Nothing else.
% ---------------------------------------------------------------------

clear all; close all  %#ok<CLALL>
set(0,'DefaultFigureWindowStyle','normal')

%% ####################################################################
%% ##  SETUP
%% ####################################################################

%% ---- [S1] Pooled input file / output ------------------------------
driveLetter = 'W:';          % drive letter for the data root and the output tree
poolFile    = fullfile(driveLetter, 'Kersten', 'Pooled', ...
                       'resonance_resultants_pool.mat');
figDir      = fullfile(driveLetter, 'Kersten', 'Figures');
figName     = 'spike_counts_exposure';
saveFig     = true;

%% ---- [S2] Inclusion mask (AND-ed; mirrors the descriptor script) --
useKS         = false;   % KS label (BC is the authoritative label gate)
useBC         = true;    % Bombcell label
useArtefact   = false;   % pulse-latency artefact gate (rate-confounded -> off)
useResponsive = true;    % ZETA pulse-responsive

%% ---- [S3] Areas (display name -> region-column value) -------------
% The pooled `region` column stores "LGN"; the thesis display label is dLGN.
% CA3 is excluded.
areaList  = { 'dLGN', "LGN" ; ...
              'V1',   "V1"  ; ...
              'CA1',  "CA1" };
areaCols  = [0.10 0.45 0.35;   % dLGN teal
             0.85 0.45 0.10;   % V1   orange
             0.20 0.30 0.55];  % CA1  navy
nArea     = size(areaList,1);

%% ---- [S4] Guard + display knobs -----------------------------------
minSpkPerBin    = 10;      % must match the descriptor and guard-sweep scripts
showPairedGuard = true;    % true : guard on Nc AND Nn (the descriptor's paired n).
                           %        Nn is the binding side under the 3 s null
                           %        window, and the one the descriptors use.
                           % false: guard on Nc only (real trials) -- hides the
                           %        effect the descriptor denominator sees.
nFreqTicks      = 8;       % how many frequency labels on the A/B index axis
greyAlpha       = 0.55;    % opacity of the below-guard overlay in [P-1 A]

fprintf('\n=== Spike-count / exposure figure (guard = %d spikes/bin, paired=%d) ===\n', ...
        minSpkPerBin, showPairedGuard);

%% ---- Load pooled file ---------------------------------------------
assert(exist(poolFile,'file')==2, 'Pooled file not found: %s', poolFile);
P = load(poolFile, 'clusterTable', 'freqList', 'Nc', 'Tc', 'Nn', 'meta');
C        = P.clusterTable;
freqList = P.freqList(:)';                 % 1 x nFreq, ascending union grid
nFreq    = numel(freqList);
fprintf('Loaded %d clusters x %d freqs (%g-%g Hz)\n', ...
        height(C), nFreq, freqList(1), freqList(end));

%% ---- Inclusion mask -----------------------------------------------
incl = true(height(C),1);
if useKS,         incl = incl & C.pass_ksLabel;     end
if useBC,         incl = incl & C.pass_bcLabel;     end
if useArtefact,   incl = incl & C.pass_notArtefact; end
if useResponsive, incl = incl & C.pass_responsive;  end
fprintf('Inclusion (KS=%d BC=%d art=%d resp=%d): %d / %d clusters.\n', ...
        useKS,useBC,useArtefact,useResponsive, sum(incl), height(C));

region = strtrim(string(C.region));

%% ####################################################################
%% ##  [1] COMPUTE : row ordering, guard mask, derived quantities
%% ####################################################################
Nc = P.Nc;  Tc = P.Tc;

% --- guard: what the downstream metric scripts actually keep ---
if showPairedGuard
    guardOK  = (Nc >= minSpkPerBin) & (P.Nn >= minSpkPerBin);
    guardLab = sprintf('N_{real} & N_{blank} \\geq %d', minSpkPerBin);
else
    guardOK  = (Nc >= minSpkPerBin);
    guardLab = sprintf('N_{real} \\geq %d', minSpkPerBin);
end
presented  = isfinite(Nc);                 % frequency existed in that recording
belowGuard = presented & ~guardOK;         % presented but dropped by the guard

% --- derived per-(cell,freq) quantities ---
rateHz     = Nc ./ max(Tc, eps);           % spikes/s during that frequency's trials
spkPerCyc  = rateHz ./ freqList;           % = Nc ./ (Tc .* f); the 1/f envelope
totSpk     = sum(Nc, 2, 'omitnan');        % total flicker spikes per cluster

% --- row ordering: group by area, sort by total spikes (descending) ---
rowIdx   = [];       % global cluster row indices, in plot order
areaSpan = zeros(nArea,2);
for a = 1:nArea
    sel = find(incl & region == areaList{a,2});
    [~, ord] = sort(totSpk(sel), 'descend');
    sel = sel(ord);
    areaSpan(a,:) = [numel(rowIdx)+1, numel(rowIdx)+numel(sel)];
    rowIdx = [rowIdx; sel]; %#ok<AGROW>
    fprintf('  %-4s : %4d cells | median total spikes %6.0f\n', ...
            areaList{a,1}, numel(sel), median(totSpk(sel),'omitnan'));
end
nCell = numel(rowIdx);

% --- reordered matrices for the heatmap ---
Ncs   = Nc(rowIdx,:);
presS = presented(rowIdx,:);
belowS= belowGuard(rowIdx,:);

% --- frequency tick positions on the index axis (log-spaced picks) ---
tickIdx = unique(round(logspace(0, log10(nFreq), nFreqTicks)));
tickLab = arrayfun(@(f) sprintf('%g', f), freqList(tickIdx), 'uni', 0);

%% ####################################################################
%% ##  [2] COMPUTE : per-area marginals
%% ####################################################################
nKept      = nan(nArea, nFreq);   % cells surviving the guard at each frequency
nAvail     = nan(nArea, nFreq);   % cells for which the frequency was presented
spcMed     = nan(nArea, nFreq);   spcLo = spcMed;  spcHi = spcMed;

for a = 1:nArea
    sel = incl & region == areaList{a,2};

    nKept(a,:)  = sum(guardOK(sel,:), 1);
    nAvail(a,:) = sum(presented(sel,:), 1);

    % Medians over ALL presented cells (NOT guard-gated): the point of panel D
    % is to show the exposure the guard is reacting to, so gating here would be
    % circular.
    S = spkPerCyc(sel,:);  S(~presented(sel,:)) = NaN;

    spcMed(a,:)  = median(S, 1, 'omitnan');
    spcLo(a,:)   = prctile(S, 25, 1);
    spcHi(a,:)   = prctile(S, 75, 1);
end

% frequency at which the MEDIAN cell drops below 1 spike/cycle (per area)
fprintf('\nMedian cell falls below 1 spike/cycle at:\n');
for a = 1:nArea
    k = find(spcMed(a,:) < 1, 1, 'first');
    if isempty(k)
        fprintf('  %-4s : never within %g-%g Hz\n', areaList{a,1}, freqList(1), freqList(end));
    else
        fprintf('  %-4s : %g Hz\n', areaList{a,1}, freqList(k));
    end
end

%% ####################################################################
%% ##  [P-1] PLOT : 4-panel exposure figure
%% ####################################################################
% Layout (tiledlayout 3x2):
%   [A] rows 1-3, col 1 : cells x frequency spike-count heatmap (full height)
%   [B] row  1,   col 2 : per-cluster total-spike histogram, per area
%   [C] row  2,   col 2 : n(f) surviving the guard      -- log-Hz axis
%   [D] row  3,   col 2 : median spikes/cycle vs freq   -- log-Hz axis
%
% Panels C and D are vertically adjacent and share a true log-Hz x-axis with
% identical ticks, so the denominator can be read directly against the exposure
% envelope. Panel A necessarily stays on a column-INDEX x-axis (imagesc cannot
% render non-uniform column spacing); its ticks carry the true frequencies.
fh = figure('Color','w','Name','Spike counts / exposure', ...
            'DefaultAxesFontSize', 9, 'DefaultTextFontSize', 9);
t  = tiledlayout(fh, 3, 2, 'TileSpacing','compact', 'Padding','loose');

% shared frequency ticks for the two log-Hz panels [C] and [D]
fTicks = [2 5 10 20 50 100 200];
fTicks = fTicks(fTicks >= freqList(1) & fTicks <= freqList(end));

% ---------------- [A] heatmap (full height, left column) -------------
axA = nexttile(t, 1, [3 1]); hold(axA,'on');
imagesc(axA, 1:nFreq, 1:nCell, log10(max(Ncs,1)), 'AlphaData', double(presS));
% grey overlay = presented but below the guard (distinct from "not presented")
greyRGB = repmat(reshape([0.55 0.55 0.55],1,1,3), nCell, nFreq);
image(axA, 1:nFreq, 1:nCell, greyRGB, 'AlphaData', greyAlpha*double(belowS));
colormap(axA, parula);
cb = colorbar(axA);
cb.Label.String = 'spikes';
cl = get(axA,'CLim');
tk = ceil(cl(1)):floor(cl(2));
set(cb, 'Ticks', tk, 'TickLabels', arrayfun(@(v) sprintf('10^{%d}',v), tk, 'uni', 0));
for a = 1:nArea                                     % area separators + labels
    if a < nArea
        yline(axA, areaSpan(a,2)+0.5, '-', 'Color',[0 0 0], 'LineWidth', 1.2);
    end
    % area label INSIDE the panel (top-left of each band) so it does not
    % collide with the colorbar sitting on the right edge
    text(axA, 1.5, areaSpan(a,1)+0.03*nCell, areaList{a,1}, ...
         'Color', areaCols(a,:), 'FontWeight','bold', ...
         'HorizontalAlignment','left', 'VerticalAlignment','top');
end
set(axA, 'YDir','reverse', 'XTick', tickIdx, 'XTickLabel', tickLab, ...
         'XLim',[0.5 nFreq+0.5], 'YLim',[0.5 nCell+0.5], 'Box','on', 'Color','w');
ylabel(axA, 'cluster (sorted by total spikes)');
title(axA, 'Spike counts (cell \times frequency)', 'FontWeight','bold');

% ---------------- [B] total-spike histogram --------------------------
axB = nexttile(t, 2); hold(axB,'on');
edges = linspace(1, ceil(log10(max(totSpk(incl)))), 40);   % log10 spikes
for a = 1:nArea
    sel = incl & region == areaList{a,2};
    v   = log10(max(totSpk(sel), 1));
    histogram(axB, v, edges, 'DisplayStyle','stairs', 'EdgeColor', areaCols(a,:), ...
              'LineWidth', 1.4, 'Normalization','count', 'DisplayName', areaList{a,1});
end
set(axB, 'Box','on');
xt = ceil(edges(1)):floor(edges(end));
set(axB, 'XTick', xt, 'XTickLabel', arrayfun(@(v) sprintf('10^{%d}',v), xt, 'uni', 0));
xlabel(axB, 'total spikes per cluster (flicker)');
ylabel(axB, 'clusters');
title(axB, 'Per-cluster yield', 'FontWeight','bold');
legend(axB, 'Location','northwest', 'Box','off');   % single key for the figure

% ---------------- [C] n(f) surviving the guard -----------------------
axC = nexttile(t, 4); hold(axC,'on');
for a = 1:nArea
    plot(axC, freqList, nAvail(a,:), ':', 'Color', areaCols(a,:), 'LineWidth', 1.0);
    plot(axC, freqList, nKept(a,:),  '-', 'Color', areaCols(a,:), 'LineWidth', 1.4);
end
set(axC, 'XScale','log', 'Box','on', 'XLim',[freqList(1) freqList(end)], ...
         'XTick', fTicks, 'XTickLabel', arrayfun(@(f) sprintf('%g',f), fTicks, 'uni', 0));
xlabel(axC, 'flicker frequency (Hz)');
ylabel(axC, 'cells');
title(axC, 'Cells entering the metric per frequency', 'FontWeight','bold');
subtitle(axC, sprintf('solid: %s   dotted: all presented', guardLab));

% ---------------- [D] spikes per cycle -------------------------------
axD = nexttile(t, 6); hold(axD,'on');
for a = 1:nArea
    plotShadedBandLU(freqList, spcMed(a,:), spcLo(a,:), spcHi(a,:), areaCols(a,:));
end
yline(axD, 1, '--', 'Color',[0.35 0.35 0.35], 'LineWidth', 1.2, ...
      'Label','1 spike / cycle', 'LabelHorizontalAlignment','left');
set(axD, 'XScale','log', 'YScale','log', 'Box','on', ...
         'XLim',[freqList(1) freqList(end)], ...
         'XTick', fTicks, 'XTickLabel', arrayfun(@(f) sprintf('%g',f), fTicks, 'uni', 0));
xlabel(axD, 'flicker frequency (Hz)');
ylabel(axD, 'spikes / cycle');
title(axD, 'Exposure per cycle = rate / f', 'FontWeight','bold');
subtitle(axD, 'the built-in 1/f envelope of |F1|   (median, IQR)');

title(t, sprintf('Raw spike counts and exposure   |   %d cells, %d mice   |   guard %s', ...
                 nCell, numel(unique(C.mouseID(incl))), guardLab), ...
      'FontWeight','bold', 'FontSize', 12);

%% ---- save figure ---------------------------------------------------
if saveFig
    if ~exist(figDir,'dir'), mkdir(figDir); end
    set(fh, 'Units','centimeters', 'Position',[2 2 40 26], 'PaperPositionMode','auto');
    drawnow;   % force label extents to settle before the tight-box export
    % exportgraphics 'Padding' is R2023b+; on older MATLAB we instead give the
    % layout a small inset so the tight bounding box does not clip labels.
    t.Padding = 'loose';
    fp = get(fh,'Position');
    set(fh, 'Position', fp + [0 0 1.2 1.2]);   % ~1 cm extra canvas each way
    drawnow;
    exportgraphics(fh, fullfile(figDir, [figName '.png']), 'Resolution', 300);
    savefig(fh, fullfile(figDir, [figName '.fig']));
    fprintf('\nFigure -> %s\n', fullfile(figDir, figName));
end

%% ####################################################################
%% ##  HELPERS  (local functions; keep at end of file)
%% ####################################################################
function plotShadedBandLU(x, mid, lo, hi, col)
    % Shaded band between explicit LOWER and UPPER bounds (e.g. quartiles),
    % with the central line on top. Sibling of plotShadedBand (mean +/- SEM);
    % separate because a median/IQR band is asymmetric and cannot be passed
    % as a single SE vector. NaN points are dropped so fill() stays a valid
    % polygon. Works on a log x-axis (set XScale AFTER calling).
    x = x(:); mid = mid(:); lo = lo(:); hi = hi(:);
    ok = isfinite(x) & isfinite(mid) & isfinite(lo) & isfinite(hi);
    x = x(ok); mid = mid(ok); lo = lo(ok); hi = hi(ok);
    if isempty(x), return; end
    fill([x; flipud(x)], [lo; flipud(hi)], col, ...
         'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(x, mid, '-', 'Color', col, 'LineWidth', 1.3);
end
