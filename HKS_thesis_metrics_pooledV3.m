%% HKS_thesis_metrics_pooledV3.m
% Thesis Figure 3 source: the three pooled main-metric curves --
%   F0    spontaneous-corrected rate (spk/s)
%   |F1|  null-subtracted first-harmonic amplitude (spk/s)
%   PPC   Vinck 2010, ungated, spike-guarded
% produced as THREE SEPARATE per-area figures (V1 / dLGN / CA1) so they can be
% stitched into one multi-panel figure by hand later.
%
% Also carries the blank-null PPC diagnostic [N], which is the evidence behind a
% decision the rest of the project depends on: PPC is count-debiased, so
% E[PPC] = 0 under a uniform null and the blank curve should sit on zero. Where
% it does, null subtraction buys no bias correction and only adds variance --
% which is why the descriptor script plots RAW PPC and uses the difference only
% for the paired test.
%
% ---------------------------------------------------------------------
% INPUT   <DRIVE>\Kersten\Pooled\resonance_resultants_pool.mat
%         No spikes, no QC files, no recSpec loop. Requires meta.dblRealDur.
% OUTPUT  <DRIVE>\Kersten\Figures\Thesis\
%           POOLED_F0_perarea_linkedV1CA1.png
%           POOLED_F1_perarea_linkedV1CA1.png
%           POOLED_PPC_perarea_linkedV1CA1.png
%           POOLED_PPCnull_perarea.png        (the [N] diagnostic)
% ---------------------------------------------------------------------
%
% ---------------------------------------------------------------------
% CURRENT BEHAVIOUR
% ---------------------------------------------------------------------
%   * V1 and CA1 SHARE a y-axis within each metric row; dLGN keeps its own.
%     Autoscaling every panel independently inflated the CA1 noise to fill its
%     frame, which argued against the claim that CA1 shows an absence. Putting
%     CA1 on the dLGN scale instead flattens V1 into an unreadable line -- dLGN
%     is roughly 5x larger in both PPC and |F1| -- so only the two panels within
%     a few-fold of one another are linked. Every panel keeps its own TICK
%     labels, which is what makes the split unambiguous.
%   * A grey blank-window noise band sits behind each curve, in one of three
%     modes per metric (see nullBandMode in [S1]): 'direct' for PPC, where a
%     matched blank estimate of the same quantity exists; 'diff' for Delta|F1|,
%     which is already a real-minus-null difference and so has no second null
%     window -- the band is centred on zero and widened by sqrt(2), a
%     conservative reference scale rather than a formal null; 'none' for F0,
%     which arrives already spontaneous-corrected and whose uncorrected blank
%     rate the pool does not carry.
%   * ONE y-label per row, on the panel named by ylabelPanel. Tick labels are
%     still drawn on every panel and already mark the dLGN / V1-CA1 scale split,
%     so repeated axis text was the larger cost. If the composite reorders the
%     panels, move ylabelPanel with it.
%   * Panel titles are "area (n = X)" only -- the metric name lives in the
%     y-label. `tag` survives as the figure window name and the filename stub.
%   * The legend sits on the panel named by legendPanel and names ONLY the grey
%     band; the coloured curve is identified by the panel title. The grey line is
%     drawn heavier than the data curve so the swatch survives reduction into the
%     composite. legendPanel = 0 drops the legend and leaves the band to the
%     caption.
%   * Blank-window PPC is computed in [1] rather than in [N], so the band on the
%     figure and the numbers in the diagnostic table cannot drift apart.
%   * F0 and |F1| are BOTH spk/s here, so the two panels are Fourier-consistent
%     and directly comparable. The per-cycle and per-period views live in the
%     supplement script.
%   * F0 divides matMuPeriod by meta.dblRealDur, NOT meta.dblStimDurPerRec.
%     matMuPeriod counts excess spikes over the ANALYSIS window (the last
%     dblRealDur s of the stimulus) -- the same window that produced Zc/Nc/Tc.
%     Dividing by the full ~4.99 s trial would under-report F0 by roughly 0.60
%     AND put F0 and |F1| on different epochs, which is exactly what these
%     panels claim they do not. An assert refuses to run on a pool that predates
%     the field.
%   * PPC is UNGATED (region mask only). A driven-cell gate was a vector-strength
%     necessity; PPC is already count-debiased and does not need one. The
%     per-bin spike guard (Nc >= minSpkPerBin) is kept, applied as a NaN mask
%     before averaging.
%   * |F1| is ungated AND unguarded (self-attenuating), and F0 carries no spike
%     guard either since it is a rate. Where the blank |F1| is undefined it is
%     treated as zero rather than masked, so an |F1| point survives where the
%     matching PPC point would be dropped.
%   * Pooling is a flat cluster mean +/- SEM across clusters. That SEM is for
%     VISUALISATION, not inference -- the animal-level bootstrap is the
%     inferential unit and lives in the descriptor script.
%   * Inclusion is frozen at reduction time (meta.qcToggles): every pooled
%     cluster already passed the pool's BC and responsiveness gates, so no
%     further inclusion mask is applied here -- region split only.
%   * Figures export via exportgraphics at exportDPI, not saveas: saveas
%     rasterises at screen resolution and the labels disappeared once panels were
%     scaled down into the composite.
%
% ---------------------------------------------------------------------
% RUN ORDER    HKS_pool_resultantsV3 -> THIS
% DEPENDS ON   Statistics Toolbox (prctile, corrcoef)
% ---------------------------------------------------------------------

clear; close all
set(0,'DefaultFigureWindowStyle','normal')

%% ####################################################################
%% ##  [S1] SETUP
%% ####################################################################
DRIVE     = 'W:';                    % drive letter for the data root and the output tree
poolFile  = fullfile(DRIVE,'Kersten','Pooled','resonance_resultants_pool.mat');
strFigDir = fullfile(DRIVE,'Kersten','Figures','Thesis');
if ~isempty(strFigDir) && ~exist(strFigDir,'dir'), mkdir(strFigDir); end
recName   = 'POOLED';                               % figure-filename stub

% ---- downstream knobs the pooler deliberately left FREE ----
minSpkPerBin   = 10;             % per-bin guard for PPC (NaN-masked before averaging)
doNullSubtract = true;           % subtract matched blank-null |F1| floor
f0Unit         = 'rate';         % 'rate' = spk/s (divide by meta.dblRealDur)
                                 % 'period' = spikes per analysis window

% ---- area set (shared across all three figures for a consistent grid) ----
areaName = {'V1','LGN','CA1'};                      % pooled region column labels
areaDisp = {'V1','dLGN','CA1'};                     % thesis display labels
areaCols = {[0 0.45 0.74],[0.85 0.33 0.10],[0.47 0.67 0.19]};
nArea    = numel(areaName);

% ---- y-axis linking -----------------------------------------------------
% Panels named here share one y-axis within each metric row. dLGN is left out
% deliberately: it is ~5x larger in PPC and ~5x larger in |F1|, and forcing V1
% onto that scale destroys the shape information the figure exists to show.
linkAreas = {'V1','CA1'};
linkIdx   = find(ismember(areaName, linkAreas));

% ---- blank-window noise band --------------------------------------------
% strNullRef, if non-empty, is appended to the legend entry as a pointer to the
% full null curves. Left EMPTY on purpose: the supplement numbering is not
% settled, and naming a specific supplement figure here would lock that layout
% in before it exists. Fill it in only once the numbering is final.
showNullBand = true;
strNullRef   = '';               % e.g. 'Supp. Fig. S1' once the supplement is fixed
%
% nullBandMode is per metric and the distinction matters for the caption:
%   'direct' -- the band IS the blank-window metric, mean +/- SEM across
%               clusters, in the metric's own units. Used for PPC, where a
%               matched blank estimate of the same quantity exists.
%   'diff'   -- the plotted quantity is already a real-minus-null difference,
%               so no second null window exists to form a null-minus-null.
%               The band is centred on zero with a width of sqrt(2) x SEM of
%               the blank term, i.e. the spread of a difference between two
%               independent draws of that size. Clusters are shared between
%               the two terms so the covariance is positive and the true
%               spread is smaller, which makes this a conservative reference
%               scale rather than a formal null distribution. Say so in the
%               caption, or set nullBandF1 = false to drop it.
%   'none'   -- F0 is already spontaneous-corrected inside the reduction and
%               the pool does not carry the uncorrected blank rate, so there
%               is nothing honest to draw here.
nullBandF1 = true;               % set false to drop the Δ|F1| band entirely

% ---- panel furniture ----------------------------------------------------
% Which panel in the 1 x nArea row carries the shared y-label, and which one
% carries the legend. Both are panel indices into areaName (1 = V1, 2 = dLGN,
% 3 = CA1); 0 disables. Tick labels stay on every panel regardless, since dLGN
% is not on the V1/CA1 scale.
ylabelPanel = 1;                 % leftmost only
legendPanel = nArea;             % CA1: flattest curve, so the most free space.
                                 % The legend now names ONLY the grey blank
                                 % band -- the coloured curve is identified by
                                 % the panel title, so it is kept out of the
                                 % legend rather than repeating "CA1" twice.

% ---- typography and export ----------------------------------------------
fs.axis   = 12;                  % tick labels
fs.label  = 13;                  % x / y labels
fs.title  = 13;                  % panel titles
fs.legend = 11;
fs.lw     = 1.6;                 % mean-curve line width
fs.lwNull = 2.0;                 % grey blank-window line: heavier than the
                                 % curve so the legend swatch stays legible
                                 % once the panel is scaled into the composite
figPos    = [100 100 1300 430];  % on-screen size; drives the exported aspect
exportDPI = 300;

%% ####################################################################
%% ##  [S2] LOAD the reduction
%% ####################################################################
assert(exist(poolFile,'file')==2, 'Pool file not found: %s', poolFile);
P = load(poolFile);
C = P.clusterTable;

Zc = P.Zc; Nc = P.Nc; Tc = P.Tc;
Zn = P.Zn; Nn = P.Nn; Tn = P.Tn;
matMuPeriod = P.matMuPeriod;
freqList    = P.freqList(:)';                        % row
nFreq       = numel(freqList);
nClust      = height(C);

region  = strtrim(string(C.region));
recKey  = strtrim(string(C.recID));

areaMask = false(nClust,nArea);
for ai = 1:nArea, areaMask(:,ai) = (region == areaName{ai}); end
fprintf('Pooled bundle: %d clusters x %d freqs from %d recording(s).\n', ...
        nClust, nFreq, numel(P.meta.recNames));
fprintf('Area masks:'); for ai=1:nArea, fprintf(' %s=%d',areaName{ai},sum(areaMask(:,ai))); end
fprintf('  (unassigned=%d)\n', sum(~any(areaMask,2)));

%% ####################################################################
%% ##  [1] COMPUTE : the three metric matrices  (nClust x nFreq)
%% ####################################################################
% ---- PPC (Vinck 2010), spike-guarded ----
VSc  = abs(Zc);
PPCc = (Nc .* VSc.^2 - 1) ./ max(Nc - 1, 1);
relMask = Nc >= minSpkPerBin;
PPCc(~relMask) = NaN;                               % per-bin guard (kept)

% ---- blank-window PPC, same estimator and same guard -------------------
% Computed here rather than in [N] because the panels in [P-1..3] draw it as
% the noise band. [N] reuses these matrices, so the band on the figure and
% the numbers in the diagnostic table can never drift apart.
VSn  = abs(Zn);
PPCn = (Nn .* VSn.^2 - 1) ./ max(Nn - 1, 1);
PPCn(Nn < minSpkPerBin) = NaN;

% strict pairing for dPPC: both sides must survive the guard, so the
% difference is never real-minus-nothing (matches the descriptor script)
bothOK   = (Nc >= minSpkPerBin) & (Nn >= minSpkPerBin);
PPCrealP = PPCc;  PPCrealP(~bothOK) = NaN;
PPCnullP = PPCn;  PPCnullP(~bothOK) = NaN;
dPPC     = PPCrealP - PPCnullP;

% ---- |F1| amplitude in spk/s, null-subtracted, unguarded ----
F1abs     = 2 .* Nc .* VSc ./ max(Tc, eps);         % spikes/s
F1absNull = 2 .* Nn .* abs(Zn) ./ max(Tn, eps);  F1absNull(isnan(F1absNull)) = 0;
if doNullSubtract, F1plot = F1abs - F1absNull; else, F1plot = F1abs; end
% (no /freqList here: 'rate' unit. per-cycle view is a supplement.)

% ---- F0 : spontaneous-corrected rate (spk/s), unguarded ----
% matMuPeriod counts excess spikes over the ANALYSIS window (the last
% dblRealDur s of the stimulus) -- the same window that produced Zc/Nc/Tc. Dividing by dblStimDurPerRec (the full ~4.99 s trial) would
% under-report F0 by dblRealDur/dblStimDur ~ 0.60 AND would put F0 and |F1| on
% different epochs, which is precisely what these panels claim they are not.
assert(isfield(P.meta,'dblRealDur'), ...
    ['Pool predates the dblRealDur field. Re-pool before plotting: the F0 ' ...
     'divisor and the Fourier window must be the same.']);
if strcmpi(f0Unit,'rate')
    F0plot = matMuPeriod ./ P.meta.dblRealDur;
    f0Ylab = 'F0  spont-corrected rate (spk/s)';
else
    F0plot = matMuPeriod;
    f0Ylab = sprintf('F0  spont-corrected spikes / %.1f s', P.meta.dblRealDur);
end

%% ####################################################################
%% ##  [S3] SHARED PLOT OPTIONS
%% ####################################################################
% Legend entry for the grey band. Bare 'blank window' unless strNullRef has
% been filled in, so the figure carries no claim about supplement numbering.
if isempty(strNullRef)
    nullName = 'blank window';
else
    nullName = sprintf('blank window (%s)', strNullRef);
end
optBase = struct( ...
    'linkIdx',     linkIdx, ...
    'ylabelPanel', ylabelPanel, ...
    'legendPanel', legendPanel, ...
    'fs',        fs, ...
    'figPos',    figPos, ...
    'exportDPI', exportDPI, ...
    'nullM',     [], ...
    'nullMode',  'none', ...
    'nullName',  nullName);

%% ####################################################################
%% ##  [P-1] F0 per area  (spk/s, ungated)
%% ####################################################################
% No null band: F0 arrives from the reduction already spontaneous-corrected
% and the uncorrected blank rate is not carried in the pool.
optF0 = optBase;
perAreaFigure(freqList, F0plot, areaMask, areaDisp, areaCols, ...
    f0Ylab, 'F0', true, ...
    fullfile(strFigDir, sprintf('%s_F0_perarea_linkedV1CA1.png', recName)), optF0);

%% ####################################################################
%% ##  [P-2] |F1| per area  (spk/s, null-subtracted, ungated)
%% ####################################################################
optF1 = optBase;
if showNullBand && nullBandF1
    optF1.nullM    = F1absNull;
    optF1.nullMode = 'diff';        % see the nullBandMode notes in [S1]
end
perAreaFigure(freqList, F1plot, areaMask, areaDisp, areaCols, ...
    '\Delta|F1| (spikes/s)', '\Delta|F1|', true, ...
    fullfile(strFigDir, sprintf('%s_F1_perarea_linkedV1CA1.png', recName)), optF1);

%% ####################################################################
%% ##  [P-3] PPC per area  (Vinck 2010, UNGATED, spike-guarded)
%% ####################################################################
optPPC = optBase;
if showNullBand
    optPPC.nullM    = PPCnullP;     % genuine matched blank estimate
    optPPC.nullMode = 'direct';
end
perAreaFigure(freqList, PPCc, areaMask, areaDisp, areaCols, ...
    'PPC (Vinck 2010)', 'PPC', true, ...    % tag = window name + filename stub
    fullfile(strFigDir, sprintf('%s_PPC_perarea_linkedV1CA1.png', recName)), optPPC);

fprintf('\nDone. Three figures written to %s\n', strFigDir);
fprintf('y-axis linked across: %s   (dLGN independent)\n', strjoin(linkAreas, ' + '));

%% ####################################################################
%% ##  [N] BLANK-NULL PPC DIAGNOSTIC
%% ####################################################################
% PPC is count-debiased, so E[PPC]=0 under a uniform null: the blank curve
% SHOULD sit on zero. Where it does not, null-subtraction is justified.
% Reports (a) blank PPC per area with animal-level bootstrap CIs, and
% (b) the empirical variance cost of subtracting it.
% PPCn / PPCnullP / dPPC are built in [1] and shared with the panels above.

nBoot = 1000;
rng(42);

% ---- [N.1] animal vector (resampling unit) ----
mouseVec = C.mouseID;
if ~isnumeric(mouseVec), mouseVec = double(string(mouseVec)); end
mouseVec = mouseVec(:);

% ---- [P-N1] blank PPC per area, bootstrap CI, real curve for scale ----
figure('Color','w','Name','Blank-null PPC per area','Position',figPos);
axN = gobjects(1,nArea);
ylN = nan(nArea,2);
for ai = 1:nArea
    rows = find(areaMask(:,ai));
    axN(ai) = subplot(1,nArea,ai); hold on
    if ~isempty(rows)
        [muN, loN, hiN] = bootAnimalCurve(PPCnullP, rows, mouseVec, nBoot);
        muR = mean(PPCrealP(rows,:), 1, 'omitnan');

        plot(freqList, muR, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2, ...
             'DisplayName','real');
        fillOK = isfinite(loN) & isfinite(hiN);
        fill([freqList(fillOK) fliplr(freqList(fillOK))], ...
             [loN(fillOK) fliplr(hiN(fillOK))], areaCols{ai}, ...
             'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
        plot(freqList, muN, '-', 'Color', areaCols{ai}, 'LineWidth', fs.lw, ...
             'DisplayName','blank null');

        % mark frequencies whose CI excludes zero
        sig = (loN > 0) | (hiN < 0);
        if any(sig)
            plot(freqList(sig), muN(sig), 'o', 'Color', areaCols{ai}, ...
                 'MarkerFaceColor', areaCols{ai}, 'MarkerSize', 4, ...
                 'HandleVisibility','off');
        end
        ylN(ai,:) = [min([loN(fillOK) muR], [], 'omitnan'), ...
                     max([hiN(fillOK) muR], [], 'omitnan')];
    end
    set(gca,'XScale','log','FontSize',fs.axis); box off
    yline(0,'k--','HandleVisibility','off');
    xlabel('Flicker frequency (Hz)','FontSize',fs.label);
    if ai == ylabelPanel                                        % [D]
        ylabel('PPC (blank null vs real)','FontSize',fs.label);
    end
    title(sprintf('%s (n = %d)', areaDisp{ai}, numel(rows)), ...  % [E]
          'FontSize',fs.title);
    if ai == legendPanel                                        % [F] was CA1
        legend('Location','best','FontSize',fs.legend,'Box','off');
    end
end
applyLinkedYLim(axN, ylN, linkIdx);
exportFigure(gcf, fullfile(strFigDir, sprintf('%s_PPCnull_perarea.png', recName)), exportDPI);

% ---- [N.2] verdict table ----
fprintf('\n=== Blank-null PPC diagnostic (guard %d spk/bin, %d boot) ===\n', ...
        minSpkPerBin, nBoot);
fprintf('%-6s %6s %10s %10s %10s %8s %8s\n', ...
        'area','nCell','med null','med real','sd ratio','r(R,N)','f CI!=0');
for ai = 1:nArea
    rows = find(areaMask(:,ai));
    if isempty(rows), continue; end
    [~, loN, hiN] = bootAnimalCurve(PPCnullP, rows, mouseVec, nBoot);
    sig = (loN > 0) | (hiN < 0);

    sdR = std(PPCrealP(rows,:), 0, 1, 'omitnan');
    sdD = std(dPPC(rows,:),     0, 1, 'omitnan');
    rRN = pairedCorrByFreq(PPCrealP(rows,:), PPCnullP(rows,:));

    fprintf('%-6s %6d %10.4f %10.4f %10.2f %8.2f %5d/%d\n', ...
        areaDisp{ai}, numel(rows), ...
        median(PPCnullP(rows,:), 'all', 'omitnan'), ...
        median(PPCrealP(rows,:), 'all', 'omitnan'), ...
        median(sdD ./ max(sdR,eps), 'omitnan'), ...
        median(rRN, 'omitnan'), sum(sig), numel(sig));
end
fprintf(['sd ratio  = median_f sd(dPPC)/sd(PPC_real). >1 means subtraction ' ...
         'adds spread.\nr(R,N)    = median across-cluster corr; positive r is ' ...
         'why the ratio can stay <1.\n']);

%% ####################################################################
%% ##  [N.3] SPONTANEOUS-WINDOW SANITY CHECK  (per-recording, optional)
%% ####################################################################
% Guarded because sQC / dblSpontWindow_s / dblStimDur belong to the
% per-recording workspace, not to this pooled script. Load a single
% recording's QC struct first if you want to run it.
if exist('sQC','var') && exist('dblSpontWindow_s','var') && exist('dblStimDur','var')
    on  = sQC.stim.flickerOn_AP(:);
    off = sQC.stim.flickerOff_AP(:);      % add pulse / noise epochs if separate
    fakeOn  = on - dblSpontWindow_s - dblStimDur;
    fakeOff = on - dblSpontWindow_s;
    nOv = sum(arrayfun(@(t) any(on < fakeOff(t) & off > fakeOn(t)), 1:numel(fakeOn)));
    fprintf('fake windows overlapping a stimulus: %d / %d\n', nOv, numel(fakeOn));
    fprintf('median onset-to-onset ITI: %.2f s  (needs > %.2f s)\n', ...
            median(diff(sort(on))), dblSpontWindow_s + dblStimDur);
else
    fprintf('\n[N.3] skipped: no per-recording sQC in the workspace.\n');
end

%% ####################################################################
%% ##  HELPERS  (local functions; keep at end of file)
%% ####################################################################

% One 1 x nArea row of panels for a single metric. `M` is nClust x nFreq;
% each panel shows the flat cluster mean +/- SEM over that area's rows.
% `zeroLine` draws a dashed y=0, for any metric whose meaningful reference is
% zero -- PPC (count-debiased, so E[PPC]=0 under a uniform null) and the
% null-subtracted |F1|.
% `opt` carries linkIdx, the typography struct, and the null-band settings.
% FORKED from the 9-argument copy in HKS_thesis_metrics_supplementV1: that one
% has no `opt` and the two are NOT interchangeable.
function perAreaFigure(freqList, M, masks, dispNames, cols, ylab, tag, zeroLine, savePath, opt)
    nArea = size(masks,2);
    fs    = opt.fs;
    figure('Color','w','Name',sprintf('%s per area', tag),'Position',opt.figPos);
    ax = gobjects(1,nArea);
    yl = nan(nArea,2);

    for ai = 1:nArea
        rows  = find(masks(:,ai));
        ax(ai) = subplot(1, nArea, ai); hold on
        extent = [];

        % ---- blank-window band first, so it sits behind the data ----
        if ~isempty(rows) && ~isempty(opt.nullM) && ~strcmpi(opt.nullMode,'none')
            [muN, seN] = meanSemAcrossClusters(opt.nullM(rows,:));
            if strcmpi(opt.nullMode,'diff')
                % difference metric: centre on zero, widen to the spread of a
                % difference between two independent draws (see the header).
                muN = zeros(size(muN));
                seN = sqrt(2) .* seN;
            end
            plotShadedBand(freqList, muN, seN, [0.55 0.55 0.55], 0.35, fs.lwNull, opt.nullName);
            extent = [extent, muN - seN, muN + seN];
        end

        % ---- the metric itself ----
        % Empty DisplayName -> plotShadedBand hides it from the legend. The
        % area is already named in the panel title, so the legend is left to
        % identify the one thing the reader cannot otherwise name.
        if ~isempty(rows)
            [mu, se] = meanSemAcrossClusters(M(rows,:));
            plotShadedBand(freqList, mu, se, cols{ai}, 0.20, fs.lw, '');
            extent = [extent, mu - se, mu + se];
        end

        set(gca,'XScale','log','FontSize',fs.axis); box off
        if zeroLine, yline(0,'k--','HandleVisibility','off'); end
        xlabel('Flicker frequency (Hz)','FontSize',fs.label);
        if ai == opt.ylabelPanel                % [D] leftmost panel only; the
            ylabel(ylab,'FontSize',fs.label);   % tick labels still mark the
        end                                     % dLGN / V1-CA1 scale split
        title(sprintf('%s (n = %d)', dispNames{ai}, numel(rows)), ...
              'FontSize',fs.title);             % [E] metric name is in ylab
        if ai == opt.legendPanel && ~isempty(opt.nullM) && ~strcmpi(opt.nullMode,'none')
            legend('Location','best','FontSize',fs.legend,'Box','off');
        end

        extent = extent(isfinite(extent));
        if ~isempty(extent), yl(ai,:) = [min(extent) max(extent)]; end
    end

    applyLinkedYLim(ax, yl, opt.linkIdx);
    if ~isempty(savePath), exportFigure(gcf, savePath, opt.exportDPI); end
end

% Give the named panels one common y-axis, padded by 5% so the curves do not
% touch the frame. Panels outside linkIdx keep whatever MATLAB chose.
function applyLinkedYLim(ax, yl, linkIdx)
    linkIdx = linkIdx(linkIdx >= 1 & linkIdx <= numel(ax));
    if numel(linkIdx) < 2, return; end
    lo = min(yl(linkIdx,1), [], 'omitnan');
    hi = max(yl(linkIdx,2), [], 'omitnan');
    if ~isfinite(lo) || ~isfinite(hi) || hi <= lo, return; end
    pad = 0.05 * (hi - lo);
    set(ax(linkIdx), 'YLim', [lo - pad, hi + pad]);
end

function [mu, se] = meanSemAcrossClusters(mat)
    mu = mean(mat, 1, 'omitnan');
    n  = sum(~isnan(mat), 1);
    se = std(mat, 0, 1, 'omitnan') ./ sqrt(max(n, 1));
end

% ---- shaded mean +/- SEM band, mean line on top --------------------
% NaN points are dropped so fill() stays a valid polygon. alpha, lw and name are
% optional. An EMPTY name keeps the mean line out of the legend -- without this,
% legend() invents a placeholder ('data1') for any line it can see.
% Shared helper: duplicated verbatim. Edit every copy or none (see README).
function plotShadedBand(x, mu, se, col, alpha, lw, name)
    if nargin < 5 || isempty(alpha), alpha = 0.2; end
    if nargin < 6 || isempty(lw),    lw    = 1.3; end
    if nargin < 7,                   name  = '';  end
    x = x(:); mu = mu(:); se = se(:);
    ok = isfinite(x) & isfinite(mu) & isfinite(se);
    x = x(ok); mu = mu(ok); se = se(ok);
    if isempty(x), return; end
    fill([x; flipud(x)], [mu - se; flipud(mu + se)], col, ...
         'FaceAlpha', alpha, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    if isempty(name)
        plot(x, mu, '-', 'Color', col, 'LineWidth', lw, 'HandleVisibility','off');
    else
        plot(x, mu, '-', 'Color', col, 'LineWidth', lw, 'DisplayName', name);
    end
end
function exportFigure(figHandle, savePath, dpi)
    try
        exportgraphics(figHandle, savePath, 'Resolution', dpi, ...
                       'BackgroundColor', 'white');
    catch
        saveas(figHandle, savePath);
    end
end

% Hierarchical bootstrap (animals, then clusters within animal) on the flat
% cluster-mean curve. Same resampling unit as the descriptor script.
function [mu, lo, hi] = bootAnimalCurve(M, rows, mouseVec, nBoot)
    mu = mean(M(rows,:), 1, 'omitnan');
    m  = mouseVec(rows);
    uM = unique(m);
    nM = numel(uM);
    idxByMouse = arrayfun(@(x) find(m == x), uM, 'uni', 0);
    B = nan(nBoot, size(M,2));
    for b = 1:nBoot
        pick = randi(nM, nM, 1);
        rr = cell(nM,1);
        for k = 1:nM
            src   = idxByMouse{pick(k)};
            rr{k} = src(randi(numel(src), numel(src), 1));
        end
        B(b,:) = mean(M(rows(vertcat(rr{:})),:), 1, 'omitnan');
    end
    lo = prctile(B, 2.5,  1);
    hi = prctile(B, 97.5, 1);
end

% Across-cluster Pearson r between real and null, computed per frequency.
function r = pairedCorrByFreq(A, Bm)
    nF = size(A,2);
    r  = nan(1,nF);
    for fi = 1:nF
        a = A(:,fi); b = Bm(:,fi);
        ok = isfinite(a) & isfinite(b);
        if sum(ok) > 3, c = corrcoef(a(ok), b(ok)); r(fi) = c(1,2); end
    end
end
