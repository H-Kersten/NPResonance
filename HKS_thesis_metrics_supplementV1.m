%% HKS_thesis_metrics_supplementV1.m
% SUPPLEMENT to Figure 2: the same F0 and |F1| pooled curves, shown in the two
% NON-RATE units (the main figure carries F0/s and |F1|/s). Four individual
% per-area figures (V1 / dLGN / CA1), for stitching by hand:
%       F0   per period   (spont-corrected spikes / window)   <- native pool unit
%       F0   per cycle    (spont-corrected spikes / cycle)
%       |F1| per cycle    (modulated spikes / cycle)           <- 1/f-tilted view
%       |F1| per period   (modulated spikes / window)
%
% WHY THESE UNITS
%   PPC is dimensionless -> one shape, no unit choice. F0 and |F1| are the only
%   rate-carrying metrics, so only they have a per-second / per-cycle /
%   per-period choice. |F1| per cycle divides the spk/s amplitude by f (a
%   monotone 1/f tilt), which turns a band-pass rate curve into an apparent
%   low-pass -- the same data, opposite qualitative shape. That unit-dependence
%   is exactly why PPC is the locked primary; this supplement is where it is
%   shown, not hidden.
%
% CONVERSIONS  (W = P.meta.dblRealDur, the ANALYSIS window -- 3.0 s, NOT the
%   ~4.99 s trial. matMuPeriod counts excess spikes over W, and Tc = nTrials*W,
%   so every unit conversion below must use W. Using dblStimDurPerRec here
%   understates F0_cycle by W/stimDur ~ 0.60 and puts the F0 and |F1|
%   "per window" panels on different epochs.):
%     F0_window = matMuPeriod                       [spikes / W]      (native)
%     F0_cycle  = matMuPeriod ./ (W .* f)           [spikes / cycle]  (= F0_rate/f)
%     F1_rate   = 2*Nc*VS/Tc  (null-subtracted)     [spikes / s]
%     F1_cycle  = F1_rate ./ f                      [mod. spikes / cycle]
%     F1_window = F1_rate .* W                      [mod. spikes / W]
%
%   Order is immaterial for the null subtraction: real and null share each
%   factor, so converting the already-subtracted amplitude equals subtracting
%   the converted ones.
%
%   Output FILENAMES keep the "perperiod" stub for continuity with the existing
%   figure-stitching workflow; the axis labels are the accurate ones.
%
% Gating matches the main script: F0 and |F1| are ungated AND unguarded (rates /
% counts, self-attenuating). Pooling is flat cluster mean +/- SEM across clusters.
%
% ---------------------------------------------------------------------
% INPUT   <DRIVE>\Kersten\Pooled\resonance_resultants_pool.mat
%         Requires meta.dblRealDur; every conversion above depends on it.
% OUTPUT  <DRIVE>\Kersten\Figures\Thesis\
%           POOLED_F0_perperiod_perarea.png
%           POOLED_F0_percycle_perarea.png
%           POOLED_F1_percycle_perarea.png
%           POOLED_F1_perperiod_perarea.png
% RUN ORDER    HKS_pool_resultantsV3 -> THIS
%              (companion to HKS_thesis_metrics_pooledV1, which carries the
%               per-second panels and the blank-null PPC diagnostic)
% ---------------------------------------------------------------------

clear; close all
set(0,'DefaultFigureWindowStyle','docked')

%% ####################################################################
%% ##  [S1] SETUP
%% ####################################################################
DRIVE     = 'W:';                    % drive letter for the data root and the output tree
poolFile  = fullfile(DRIVE,'Kersten','Pooled','resonance_resultants_pool.mat');
strFigDir = fullfile(DRIVE,'Kersten','Figures','Thesis');
if ~isempty(strFigDir) && ~exist(strFigDir,'dir'), mkdir(strFigDir); end
recName   = 'POOLED';

doNullSubtract = true;           % subtract matched blank-null |F1| floor

areaName = {'V1','LGN','CA1'};
areaDisp = {'V1','dLGN','CA1'};
areaCols = {[0 0.45 0.74],[0.85 0.33 0.10],[0.47 0.67 0.19]};
nArea    = numel(areaName);

% which supplement panels to emit (set any to false to skip)
do_F0_period = true;
do_F0_cycle  = true;
do_F1_cycle  = true;
do_F1_period = true;

%% ####################################################################
%% ##  [S2] LOAD the reduction
%% ####################################################################
assert(exist(poolFile,'file')==2, 'Pool file not found: %s', poolFile);
P = load(poolFile);
C = P.clusterTable;

Zc = P.Zc; Nc = P.Nc; Tc = P.Tc;
Zn = P.Zn; Nn = P.Nn; Tn = P.Tn;
matMuPeriod = P.matMuPeriod;
freqList    = P.freqList(:)';                        % 1 x nFreq (row)
nFreq       = numel(freqList);
nClust      = height(C);

region  = strtrim(string(C.region));
recKey  = strtrim(string(C.recID));

areaMask = false(nClust,nArea);
for ai = 1:nArea, areaMask(:,ai) = (region == areaName{ai}); end

% Analysis window (scalar, identical across recordings). The window is a fixed
% span inside the trial, not the whole trial, so this is a scalar rather than a
% per-cluster join against dblStimDurPerRec.
assert(isfield(P.meta,'dblRealDur'), ...
    ['Pool predates the dblRealDur field. Re-pool before plotting: every unit ' ...
     'conversion here depends on the analysis window.']);
W = P.meta.dblRealDur;
fprintf('Analysis window W = %.2f s (trial length %.2f s).\n', ...
        W, max(P.meta.dblStimDurPerRec));

fprintf('Pooled bundle: %d clusters x %d freqs from %d recording(s).\n', ...
        nClust, nFreq, numel(P.meta.recNames));

%% ####################################################################
%% ##  [1] COMPUTE : native quantities + unit variants
%% ####################################################################
% ---- |F1| amplitude in spk/s, null-subtracted (the native rate quantity) ----
VSc       = abs(Zc);
F1abs     = 2 .* Nc .* VSc ./ max(Tc, eps);
F1absNull = 2 .* Nn .* abs(Zn) ./ max(Tn, eps);  F1absNull(isnan(F1absNull)) = 0;
if doNullSubtract, F1_rate = F1abs - F1absNull; else, F1_rate = F1abs; end

% ---- unit variants (W scalar, freqList is 1 x nFreq) ----
F0_period = matMuPeriod;                              % spikes / W (native)
F0_cycle  = matMuPeriod ./ (W .* freqList);           % spikes / cycle
F1_cycle  = F1_rate ./ freqList;                      % modulated spikes / cycle
F1_period = F1_rate .* W;                             % modulated spikes / W

%% ####################################################################
%% ##  [P] PLOTS : one 1 x nArea figure per (metric, unit)
%% ####################################################################
if do_F0_period
    perAreaFigure(freqList, F0_period, areaMask, areaDisp, areaCols, ...
        sprintf('F0  spont-corr. spikes / %.1f s window', W), 'F0/window', false, ...
        fullfile(strFigDir, sprintf('%s_F0_perperiod_perarea.png', recName)));
end
if do_F0_cycle
    perAreaFigure(freqList, F0_cycle, areaMask, areaDisp, areaCols, ...
        'F0  spont-corr. spikes / cycle', 'F0/cycle', false, ...
        fullfile(strFigDir, sprintf('%s_F0_percycle_perarea.png', recName)));
end
if do_F1_cycle
    perAreaFigure(freqList, F1_cycle, areaMask, areaDisp, areaCols, ...
        '\Delta|F1|  modulated spikes / cycle', '\Delta|F1|/cycle', true, ...
        fullfile(strFigDir, sprintf('%s_F1_percycle_perarea.png', recName)));
end
if do_F1_period
    perAreaFigure(freqList, F1_period, areaMask, areaDisp, areaCols, ...
        sprintf('\\Delta|F1|  modulated spikes / %.1f s window', W), '\Delta|F1|/window', true, ...
        fullfile(strFigDir, sprintf('%s_F1_perperiod_perarea.png', recName)));
end

fprintf('\nDone. Supplement figures written to %s\n', strFigDir);

%% ####################################################################
%% ##  HELPERS  (shared with HKS_thesis_metrics_pooledV1.m; keep in step)
%% ####################################################################

% One 1 x nArea row of panels for a single metric. `M` is nClust x nFreq;
% each panel shows the flat cluster mean +/- SEM over that area's rows.
% `zeroLine` draws a dashed y=0, for any metric whose meaningful reference is
% zero -- PPC (count-debiased, so E[PPC]=0 under a uniform null) and the
% null-subtracted |F1|.
% Shared helper, with ONE known divergence: the pooled-metrics copy carries a
% disabled linkaxes block. Keep the two in step (see README).
function perAreaFigure(freqList, M, masks, disp, cols, ylab, tag, zeroLine, savePath)
    nArea = size(masks,2);
    figure('Color','w','Name',sprintf('%s per area', tag));
    for ai = 1:nArea
        rows = find(masks(:,ai));
        subplot(1, nArea, ai); hold on
        if ~isempty(rows)
            [mu, se] = meanSemAcrossClusters(M(rows,:));
            plotShadedBand(freqList, mu, se, cols{ai});
        end
        set(gca,'XScale','log'); box off
        if zeroLine, yline(0,'k--','HandleVisibility','off'); end
        xlabel('Flicker frequency (Hz)');
        if ai == 1, ylabel(ylab); end
        title(sprintf('%s  %s  (n=%d)', disp{ai}, tag, numel(rows)));
    end
    if ~isempty(savePath), saveas(gcf, savePath); end
end

% ---- mean +/- SEM across clusters ----------------------------------
% Shared helper: duplicated verbatim. Edit every copy or none (see README).
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
