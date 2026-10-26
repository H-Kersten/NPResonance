%% HKS_thesis_descriptors_PPCV2.m
% Thesis Figure 4 source: PPC following-descriptors (f_max / f_half / f_cutoff)
% drawn at THREE grouping levels, each as its OWN single-axes figure:
%     'area'        -> dLGN / V1 / CA1
%     'v1layer'     -> V1 laminar split (L2/3 / L4 / L5 / L6)
%     'v1celltype'  -> V1 NS / BS (narrow / broad-spiking)
% Mode 4 (v1layertype) is skipped: layer x cell-type is too thin for
% animal-level CIs.
%
% Stage 4: the last script on the critical path.
%
% ---------------------------------------------------------------------
% INPUTS
% ---------------------------------------------------------------------
%   <DRIVE>\Kersten\Pooled\resonance_resultants_pool.mat
%       Zc/Nc and Zn/Nn drive the metrics; clusterTable supplies region,
%       v1Layer, cellType, mouseID, recID and the pass_* columns. No spikes.
%       v1Layer and cellType are read from the pooled columns when present and
%       derived in place otherwise, so the pool stays usable if either is missing.
%
% ---------------------------------------------------------------------
% OUTPUT  ->  <DRIVE>\Kersten\Figures\Thesis\
%                 following_descriptors_PPC_<mode>_<bootUnit>.png   (one per mode)
% ---------------------------------------------------------------------
%   allTables.<mode>   descriptor table per grouping level: group, nCell, f_max,
%                      f_half, f_cutoff (Hz) plus bootstrap CI columns.
%                      Left in the workspace; nothing is written to disk.
%   Save dir is FLAT -- the filename carries the mode, so there are no per-level
%   sub-folders here.
%
% ---------------------------------------------------------------------
% CURRENT BEHAVIOUR
% ---------------------------------------------------------------------
%   * TWO curves, guarded differently on purpose. Mcurve is RAW PPC guarded on
%     the real side only, and drives the plotted curve and the geometric
%     descriptors f_max / f_half. Dprim is the paired real-minus-blank
%     difference guarded on BOTH sides, and drives only the signrank test behind
%     f_cutoff. Blank PPC is statistically zero, so subtracting it buys no bias
%     correction and only adds spread, worst where the signal is weakest -- but a
%     per-cluster empirical null is still more defensible than testing against an
%     analytic zero, which is why the paired test keeps it.
%   * NULL GATE. If no frequency survives the paired test, all three descriptors
%     report NaN (n.d.), never a floor value: max() and the half-crossing of
%     noise are meaningless, and a bare 2.0 Hz would read as real low-frequency
%     following. CA1 is the case this exists for.
%   * The bootstrap CIs are gated on the SAME rule, so a gated group returns no
%     interval. Otherwise any contrast against it compares with an undefined
%     quantity.
%   * The per-(cell,freq) spike guard is variance control, not bias control. It
%     gates real AND blank, so groups that differ in firing rate -- as NS and BS
%     do -- have their high-frequency curve computed over a shrinking,
%     rate-biased subset. A group f_cutoff gap that tracks the survivor counts
%     rather than the curve is a selection artefact; [2b] prints both counts so
%     that is visible.
%   * PPC only. The |F1| companion panel is not drawn here: PPC is the locked
%     primary, and F1 at every grouping level is more than the argument needs.
%   * Mode-independent work (load, inclusion, covariates, the per-(cell,freq)
%     PPC difference) is hoisted above the mode loop and runs once. rng is
%     reseeded per mode, so each mode's bootstrap is reproducible independently
%     of how many modes ran before it.
%   * Titles are drawn only for the modes named in titleModes. The marker key
%     is printed to the console once per run as a ready-made caption line, so it
%     is stated once in the caption rather than three times inside the panels.
%   * Legend entries carry n cells and n mice per group, and a group whose
%     descriptors were gated to NaN is tagged "n.d.". Animal coverage is uneven
%     across layers, and this puts that on the face of the figure rather than
%     only in the table.
%   * Type is sized for panels that get reduced into a stitched composite, so it
%     looks oversized on screen on purpose. Export goes through exportgraphics
%     at exportDPI; saveas rasterises at screen resolution and the labels
%     disappeared once panels were scaled down.
%   * xLimHz is off by default. The flat stretch out to 243 Hz is itself a
%     result, so trimming it is a deliberate choice rather than a default.
%
% ---------------------------------------------------------------------
% RUN ORDER    HKS_pool_resultantsV3 -> THIS
% DEPENDS ON   Statistics Toolbox (signrank, prctile). No spikes, no zetatest.
% CONVENTIONS  fullfile only, no cd; SETUP -> COMPUTE -> PLOTS
% RELATED      HKS_following_descriptorsV6.m is superseded and produces different
%              numbers; its header lists the differences.
% ---------------------------------------------------------------------

clear; close all
set(0,'DefaultFigureWindowStyle','normal')

%% ####################################################################
%% ##  SETUP  (mode-independent)
%% ####################################################################

%% ---- [S1] Pooled input file + save dir ----------------------------
driveLetter = 'W:';                                  % drive letter for the data root
poolFile    = fullfile(driveLetter, 'Kersten', 'Pooled', ...
                       'resonance_resultants_pool.mat');
strFigDir   = fullfile(driveLetter, 'Kersten', 'Figures', 'Thesis');
if ~exist(strFigDir, 'dir'), mkdir(strFigDir); end

%% ---- [S2] Modes to run (mode 4 skipped) --------------------------
modeList = {'area','v1layer','v1celltype'};

%% ---- [S3] Inclusion mask (re-applied downstream) -----------------
useKS         = false;
useBC         = true;
useArtefact   = false;
useResponsive = true;

%% ---- [S4] Layer / cell-type group definitions --------------------
layerDefs = { 'L2/3', {"L2/3"} ; ...
              'L4',   {"L4"}   ; ...
              'L5',   {"L5"}   ; ...
              'L6',   {"L6a","L6b"} };
typeDefs   = { 'NS', {"NS"} ; ...
               'BS', {"BS"} };
NS_CUTOFF_US = 400;          % peak-trough; < = NS, >= = BS (fallback derivation only)

%% ---- [S5] Metric + curve parameters ------------------------------
metricPrimary  = 'PPC';      % unused; the figure name hardcodes PPC
minSpkPerBin   = 10;         % per-(cell,freq) guard on PPC; variance control, not bias
halfMode       = 'amplitude';
switch lower(halfMode)
    case 'amplitude', halfFrac = 0.5;
    case 'power3db',  halfFrac = 1/sqrt(2);
    otherwise, error('halfMode must be amplitude | power3dB');
end

%% ---- [S6] f_cutoff test parameters -------------------------------
cutoffAlpha    = 0.05;       % BH-FDR level across the frequency grid
nLostRun       = 2;          % consecutive n.s. bins that define "staying lost"
minCellsForTest= 5;          % min finite paired cells to attempt signrank at a freq
gateOnSignificance = true;   % no frequency survives the paired test -> ALL three
                             % descriptors report NaN (n.d.), never a floor value

%% ---- [S7] Bootstrap parameters -----------------------------------
doBootstrap    = true;
bootUnit       = 'animal';   % 'animal' (honest/wide) | 'recording'
nBoot          = 1000;
bootCI         = [2.5 97.5]; % percentile CI
cutoffBootMethod = 'signrank';  % 'signrank' (faithful, slow) | 'curveFloor' (fast)

%% ---- [S9] Figure appearance --------------------------------------
% Sizes are set for panels that get reduced into a stitched composite, so they
% look oversized on screen on purpose.
fs.axis    = 12;             % tick labels
fs.label   = 13;             % x / y labels
fs.title   = 12;             % panel title (only drawn for titleModes)
fs.legend  = 13;             % legend text
fs.guide   = 11;             % the f_half tick label
fs.lw      = 1.6;            % mean-curve line width
figPos     = [100 100 1400 620];
exportDPI  = 300;

titleModes    = {'area'};    % modes that keep an in-panel title; {} for none
legendShowN   = true;        % append "(n=.., k mice)" and an n.d. tag
legendLoc     = 'northeast';
xLimHz        = [];          % e.g. [2 120] to trim the flat tail; [] = full grid

% Caption-ready key, printed once instead of repeated in every title.
keyLine = ['Curve: raw PPC, mean +/- SEM across cells. Marker: f_max. ' ...
           'Bracket: f_half. Dashed vertical: f_cutoff, from the paired ' ...
           'test against the blank window.'];

fprintf('\n=== Thesis PPC descriptors : modes = %s ===\n', strjoin(modeList,', '));

%% ---- Load pooled file --------------------------------------------
assert(exist(poolFile,'file')==2, 'Pooled file not found: %s', poolFile);
P = load(poolFile);
C = P.clusterTable;
freqList = P.freqList(:)';
nFreq    = numel(freqList);
fprintf('Loaded %d clusters x %d freqs from %s\n', height(C), nFreq, poolFile);

%% ---- Inclusion mask (pass_* toggles, AND-ed) ---------------------
incl = true(height(C),1);
if useKS,         incl = incl & C.pass_ksLabel;     end
if useBC,         incl = incl & C.pass_bcLabel;     end
if useArtefact,   incl = incl & C.pass_notArtefact; end
if useResponsive, incl = incl & C.pass_responsive;  end
fprintf('Inclusion (KS=%d BC=%d art=%d resp=%d): %d / %d clusters.\n', ...
        useKS,useBC,useArtefact,useResponsive, sum(incl), height(C));

region  = strtrim(string(C.region));
mouseID = C.mouseID;                     % numeric
recID   = strtrim(string(C.recID));      % string

%% ---- v1Layer covariate (pooled column, else derive from regionAllen) ----
if any(strcmp(C.Properties.VariableNames, 'v1Layer'))
    v1Layer = strtrim(string(C.v1Layer));
    fprintf('v1Layer: using pooled column.\n');
else
    assert(any(strcmp(C.Properties.VariableNames, 'regionAllen')), ...
        ['Pooled table has neither v1Layer nor regionAllen -- re-pool with ' ...
         'HKS_cluster_QC_unifiedV7 to get laminar labels.']);
    v1Layer = repmat("", height(C), 1);
    for ii = 1:height(C)
        t = regexp(string(C.regionAllen(ii)), 'layer\s+([\w/]+)', ...
                   'tokens', 'once', 'ignorecase');
        if ~isempty(t), v1Layer(ii) = "L" + string(t{1}); end
    end
    fprintf('v1Layer: derived in place from regionAllen (same parse as the QC script).\n');
end

%% ---- cellType covariate (pooled column, else derive from wfDuration_us) --
% Always needed: v1celltype is in the mode loop.
cellType = repmat("", height(C), 1);
if any(strcmp(C.Properties.VariableNames, 'cellType'))
    cellType = strtrim(string(C.cellType));
    fprintf('cellType: using pooled column.\n');
elseif any(strcmp(C.Properties.VariableNames, 'wfDuration_us'))
    wf = double(C.wfDuration_us);
    cellType(wf <  NS_CUTOFF_US) = "NS";
    cellType(wf >= NS_CUTOFF_US) = "BS";     % NaN width stays "" -> excluded
    fprintf('cellType: derived in place from wfDuration_us at %g us.\n', NS_CUTOFF_US);
else
    error(['Pooled table has neither cellType nor wfDuration_us. Rerun: ' ...
           'HKS_bombcell_batchV3 (wfDuration_us -> TSV) -> ' ...
           'HKS_cluster_QC_unified (cfg.NS_CUTOFF_US) -> HKS_pool_resultants.']);
end

%% ####################################################################
%% ##  [1] COMPUTE : per-(cluster,freq) PPC paired difference D_PPC
%% ####################################################################
% Metric formulas verbatim from HKS_fourierV4.m D.1-D.3. Mode-independent, so
% computed ONCE. (|F1| companion dropped: PPC-only script.)
Zc = P.Zc; Nc = P.Nc; Tc = P.Tc;
Zn = P.Zn; Nn = P.Nn; Tn = P.Tn;

VS_r  = abs(Zc);
PPC_r = (Nc .* VS_r.^2 - 1) ./ max(Nc - 1, 1);           % Vinck 2010
VS_n  = abs(Zn);
PPC_n = (Nn .* VS_n.^2 - 1) ./ max(Nn - 1, 1);

% per-bin spike guard (require BOTH real & null valid)
relReal = Nc >= minSpkPerBin;
relNull = Nn >= minSpkPerBin;
pairOK  = relReal & relNull;

% TWO curves, deliberately guarded differently (see the null-subtraction
% audit, 2026-07). The blank PPC is statistically zero, so subtracting it
% buys no bias correction and only adds spread -- worst where the signal is
% weakest. So:
%   Mcurve : RAW PPC, guarded on the REAL side only. Drives the plotted curve
%            and the geometric descriptors f_max / f_half (and kMax).
%   Dprim  : the paired difference, guarded on BOTH sides. Drives the paired
%            signrank test behind f_cutoff, where a per-cluster empirical null
%            is still more defensible than testing against an analytic zero.
D_PPC  = PPC_r - PPC_n;   D_PPC(~pairOK) = NaN;      % paired test -> f_cutoff
Dprim  = D_PPC;
Mcurve = PPC_r;           Mcurve(~relReal) = NaN;    % curve -> f_max / f_half
primLab = 'PPC (Vinck 2010)';
fprintf('Real-valid   (cell,freq) entries: %d / %d (guard Nc >= %d)\n', ...
        sum(relReal(:)), numel(relReal), minSpkPerBin);
fprintf('Paired-valid (cell,freq) entries: %d / %d (guard Nc,Nn >= %d)\n', ...
        sum(pairOK(:)), numel(pairOK), minSpkPerBin);

%% ####################################################################
%% ##  MODE LOOP
%% ####################################################################
allTables = struct();
for mi = 1:numel(modeList)
    groupMode = modeList{mi};
    rng(7, 'twister');                       % per-mode reproducible resamples
    fprintf('\n================ MODE: %s ================\n', groupMode);

    %% ---- [S8] Build the group set for this mode -------------------
    switch lower(groupMode)
        case 'area'
            % display label 'dLGN' maps to the stored region string "LGN"
            areaList   = { 'dLGN', 'LGN' ; 'V1', 'V1' ; 'CA1', 'CA1' };
            areaColors = {[0.85 0.33 0.10], [0 0.45 0.74], [0.47 0.67 0.19]};
            groupSel   = @(a) strcmp(region, areaList{a,2});
            groupLab   = 'area';
        case 'v1layer'
            areaList   = layerDefs(:,1);
            areaColors = {[0.20 0.55 0.90], [0.10 0.70 0.55], ...
                          [0.95 0.60 0.15], [0.60 0.30 0.65]};
            groupSel   = @(a) (region == "V1") & ismember(v1Layer, string(layerDefs{a,2}));
            groupLab   = 'V1 layer';
        case 'v1celltype'
            areaList   = typeDefs(:,1);
            areaColors = {[0.85 0.33 0.10], [0.20 0.30 0.55]};   % NS warm, BS cool
            groupSel   = @(a) (region == "V1") & ismember(cellType, string(typeDefs{a,2}));
            groupLab   = 'V1 cell type';
        otherwise
            error('modeList entries must be area | v1layer | v1celltype');
    end
    nArea = size(areaList,1);
    assert(numel(areaColors) >= nArea, 'areaColors: need >= %d colors.', nArea);

    % sanity prints (catch a stale pool / failed parse immediately)
    if strcmpi(groupMode,'v1layer')
        uL = unique(v1Layer(region == "V1" & v1Layer ~= ""));
        fprintf('V1 layers present:');
        for k = 1:numel(uL), fprintf(' %s=%d', uL(k), sum(region=="V1" & v1Layer==uL(k))); end
        fprintf('  | unlabelled V1=%d\n', sum(region=="V1" & v1Layer==""));
        assert(~isempty(uL), 'No V1 laminar labels found -- check regionAllen in the pool.');
    elseif strcmpi(groupMode,'v1celltype')
        inclV1 = incl & (region == "V1");
        fprintf('V1 cell types (included): NS=%d  BS=%d  | unlabelled V1=%d\n', ...
            sum(inclV1 & cellType=="NS"), sum(inclV1 & cellType=="BS"), ...
            sum(inclV1 & cellType==""));
    end

    %% ---- [2] Point-estimate descriptors --------------------------
    areaRows = cell(nArea,1);
    for a = 1:nArea, areaRows{a} = find(incl & groupSel(a)); end

    pt = struct('area',{},'f_max',{},'f_half',{},'f_cutoff',{}, ...
                'peak',{},'nCell',{},'halfCensored',{},'cutCensored',{});
    DareaPPC = nan(nArea, nFreq);   semPPC = nan(nArea, nFreq);

    fprintf('\n--- Point-estimate descriptors ---\n');
    for a = 1:nArea
        rows = areaRows{a};
        [muP, seP] = meanSemAcrossClusters(Mcurve(rows,:));
        DareaPPC(a,:) = muP;  semPPC(a,:) = seP;

        cd = extractCurveDescriptors(freqList, muP, halfFrac);
        [f_cut, sigVec, ~, ~, cutCens] = cutoffCrossing( ...
            freqList, Dprim(rows,:), cd.kMax, cutoffAlpha, nLostRun, minCellsForTest);

        pt(a).area         = areaList{a,1};
        pt(a).f_max        = cd.f_max;
        pt(a).f_half       = cd.f_half;
        pt(a).f_cutoff     = f_cut;
        pt(a).peak         = cd.peak;
        pt(a).nCell        = numel(rows);
        pt(a).halfCensored = cd.halfCensored;
        pt(a).cutCensored  = cutCens;
        pt(a).sigVec       = sigVec;   %#ok<*STRNU>

        % NULL GATE: if no frequency survives the paired test, the curve is
        % noise and max()/half-crossing of noise are meaningless. Report n.d.,
        % never a floor value (a bare 2.0 Hz reads as real low-freq following).
        pt(a).sigAny = any(sigVec);
        if gateOnSignificance && ~pt(a).sigAny
            pt(a).f_max = NaN;  pt(a).f_half = NaN;  pt(a).f_cutoff = NaN;
        end

        fprintf('  %-6s n=%4d | f_max=%6.1f  f_half=%6.1f%s  f_cutoff=%6.1f%s Hz%s\n', ...
            pt(a).area, pt(a).nCell, pt(a).f_max, pt(a).f_half, ...
            tern(cd.halfCensored,'*',' '), pt(a).f_cutoff, tern(cutCens,'*',' '), ...
            tern(~pt(a).sigAny,'   <- gated (n.d.)',''));
    end
    fprintf('  (* = right-censored: not reached within the tested %g-%g Hz grid)\n', ...
            min(freqList), max(freqList));

    %% ---- [2b] Guard-survival audit (selection-bias check) --------
    % minSpkPerBin gates real AND blank per (cell,freq). If groups differ in
    % firing rate -- as NS vs BS do -- the high-frequency curve is computed over
    % a shrinking, rate-biased subset of the slower group. A group f_cutoff gap
    % that tracks the survivor counts (not the curve) is a selection artefact.
    nValid     = nan(nArea, nFreq);   % paired -> gates the f_cutoff test
    nValidCurv = nan(nArea, nFreq);   % real-only -> gates the plotted curve
    for a = 1:nArea
        nValid(a,:)     = sum(isfinite(Dprim(areaRows{a},:)),  1);
        nValidCurv(a,:) = sum(isfinite(Mcurve(areaRows{a},:)), 1);
    end
    fprintf(['Guard survival at %g Hz: curve %s | paired %s (cells per group)\n'], ...
            freqList(end), mat2str(nValidCurv(:,end)'), mat2str(nValid(:,end)'));
    thinBin = any(nValid < minCellsForTest, 1);
    if any(thinBin)
        fprintf(['NOTE: %d freq bin(s) drop below minCellsForTest=%d for >=1 group; ' ...
                 'descriptors there rest on <%d cells.\n'], ...
                 sum(thinBin), minCellsForTest, minCellsForTest);
    end
    if strcmpi(groupMode,'v1celltype')
        fprintf(['NOTE (NS/BS): the minSpkPerBin sweep is COMPLETE -- see\n' ...
                 '  HKS_guard_sweep_tablesV1 for the tested values and results.\n' ...
                 '  If you change minSpkPerBin, re-check that any group f_cutoff\n' ...
                 '  gap tracks the CURVE and not nValid.\n']);
    end

    %% ---- [3] Hierarchical bootstrap -> CIs + between-group contrasts
    bootTab = struct(); contrastTab = struct();
    if doBootstrap
        switch lower(bootUnit)
            case 'animal',    unitKey = mouseID;  unitList = unique(mouseID(incl));
            case 'recording', unitKey = recID;    unitList = unique(recID(incl));
            otherwise, error('bootUnit must be animal | recording');
        end
        nUnit = numel(unitList);
        fprintf('\n--- Hierarchical bootstrap (%s, n=%d units, %d iters) ---\n', ...
                bootUnit, nUnit, nBoot);
        if nUnit < 4, warning('Only %d bootstrap units -- CIs will be very coarse.', nUnit); end

        unitRowsPerArea = cell(nArea, nUnit);
        for a = 1:nArea
            for u = 1:nUnit
                unitRowsPerArea{a,u} = areaRows{a}(unitKey(areaRows{a}) == unitList(u));
            end
        end

        floorLevel = nan(nArea,1);
        for a = 1:nArea
            if isfinite(pt(a).f_cutoff)
                [~,kc] = min(abs(freqList - pt(a).f_cutoff));
                floorLevel(a) = DareaPPC(a,kc);
            end
        end

        bMax = nan(nBoot,nArea); bHalf = nan(nBoot,nArea); bCut = nan(nBoot,nArea);
        for b = 1:nBoot
            draw = unitList(randi(nUnit, nUnit, 1));
            for a = 1:nArea
                sel = cell(numel(draw),1);
                for d = 1:numel(draw)
                    uu = find(unitList == draw(d), 1);
                    sel{d} = unitRowsPerArea{a,uu};
                end
                rows = vertcat(sel{:});
                if isempty(rows), continue; end
                mu = mean(Mcurve(rows,:), 1, 'omitnan');
                cd = extractCurveDescriptors(freqList, mu, halfFrac);
                bMax(b,a)  = cd.f_max;
                bHalf(b,a) = cd.f_half;
                switch lower(cutoffBootMethod)
                    case 'signrank'
                        fc = cutoffCrossing(freqList, Dprim(rows,:), cd.kMax, ...
                                            cutoffAlpha, nLostRun, minCellsForTest);
                    case 'curvefloor'
                        fc = crossingByFloor(freqList, mu, cd.kMax, floorLevel(a), nLostRun);
                    otherwise, error('cutoffBootMethod must be signrank | curveFloor');
                end
                bCut(b,a) = fc;
            end
        end

        % Gate the CIs on the SAME rule as the point estimates, or a gated
        % group still returns an interval and any contrast against it is a
        % comparison with an undefined quantity.
        if gateOnSignificance
            for a = 1:nArea
                if ~pt(a).sigAny
                    bMax(:,a) = NaN;  bHalf(:,a) = NaN;  bCut(:,a) = NaN;
                end
            end
        end

        descNames = {'f_max','f_half','f_cutoff'};
        descBoot  = {bMax, bHalf, bCut};
        fprintf('\n  Descriptor  ');
        for a = 1:nArea, fprintf('%-22s', areaList{a,1}); end
        fprintf('\n');
        for k = 1:3
            M = descBoot{k};
            lo = prctile(M, bootCI(1), 1); hi = prctile(M, bootCI(2), 1);
            md = median(M, 1, 'omitnan');
            fprintf('  %-8s ', descNames{k});
            for a = 1:nArea, fprintf('%6.1f [%5.1f-%5.1f] ', md(a), lo(a), hi(a)); end
            fprintf('\n');
            bootTab.(descNames{k}) = struct('median',md,'lo',lo,'hi',hi,'raw',M);
        end

        pairs = nchoosek(1:nArea, 2);
        fprintf('\n  Between-group contrasts (median diff [%.0f-%.0f%% CI]; * excludes 0):\n', ...
                bootCI(1), bootCI(2));
        for k = 1:3
            M = descBoot{k};
            for pI = 1:size(pairs,1)
                i = pairs(pI,1); j = pairs(pI,2);
                dd = M(:,i) - M(:,j);
                lo = prctile(dd, bootCI(1)); hi = prctile(dd, bootCI(2));
                md = median(dd,'omitnan');
                star = tern(~(lo<=0 && hi>=0),'*',' ');
                fprintf('   %-8s %-6s - %-6s : %6.1f [%6.1f, %6.1f] %s\n', ...
                    descNames{k}, areaList{i,1}, areaList{j,1}, md, lo, hi, star);
                fname = matlab.lang.makeValidName( ...
                    sprintf('%s_%s_%s', descNames{k}, areaList{i,1}, areaList{j,1}));
                contrastTab.(fname) = struct('median',md,'lo',lo,'hi',hi, ...
                    'excludesZero',~(lo<=0 && hi>=0), ...
                    'groupA',areaList{i,1},'groupB',areaList{j,1});
            end
        end
    end

    %% ---- assemble descriptor table for this mode -----------------
    descriptorTable = table( ...
        string({pt.area})', [pt.nCell]', [pt.f_max]', [pt.f_half]', [pt.f_cutoff]', ...
        'VariableNames', {'group','nCell','f_max','f_half','f_cutoff'});
    if doBootstrap
        descriptorTable.f_max_CI   = [bootTab.f_max.lo(:),    bootTab.f_max.hi(:)];
        descriptorTable.f_half_CI  = [bootTab.f_half.lo(:),   bootTab.f_half.hi(:)];
        descriptorTable.f_cut_CI   = [bootTab.f_cutoff.lo(:), bootTab.f_cutoff.hi(:)];
    end
    disp(descriptorTable);
    allTables.(matlab.lang.makeValidName(groupMode)) = descriptorTable;

    %% ---- [P] PPC-only figure (single axes) -----------------------
    figure('Color','w','Name',sprintf('PPC following descriptors - %s', groupLab), ...
           'Position', figPos);
    hold on
    h = gobjects(nArea,1);
    legLab = cell(nArea,1);
    for a = 1:nArea
        plotShadedBand(freqList, DareaPPC(a,:), semPPC(a,:), areaColors{a}, [], fs.lw, '');
        h(a) = plot(NaN,NaN,'-','Color',areaColors{a},'LineWidth',fs.lw);  % legend proxy
        guideCol = 0.7 * areaColors{a};

        % legend label: group, cells, animals, and whether it yielded descriptors
        legLab{a} = areaList{a,1};
        if legendShowN
            nMice = numel(unique(mouseID(areaRows{a})));
            legLab{a} = sprintf('%s (n=%d, %d mice)', legLab{a}, pt(a).nCell, nMice);
            if ~pt(a).sigAny, legLab{a} = [legLab{a} '  n.d.']; end
        end

        % f_max marker
        if isfinite(pt(a).f_max)
            [~,km] = min(abs(freqList-pt(a).f_max));
            plot(pt(a).f_max, DareaPPC(a,km), 'v', 'Color', areaColors{a}, ...
                 'MarkerFaceColor', areaColors{a}, 'MarkerSize', 6, 'HandleVisibility','off');
        end
        % f_half : L-bracket at the half-max level
        if isfinite(pt(a).f_half) && isfinite(pt(a).peak) && pt(a).peak > 0
            yHalf = halfFrac * pt(a).peak;
            plot([pt(a).f_max pt(a).f_half], [yHalf yHalf], ':', ...
                 'Color', guideCol, 'LineWidth', 1.6, 'HandleVisibility','off');
            plot([pt(a).f_half pt(a).f_half], [0 yHalf], ':', ...
                 'Color', guideCol, 'LineWidth', 1.6, 'HandleVisibility','off');
            text(pt(a).f_half, 0, '  f_{half}', 'Color', guideCol, 'FontSize', fs.guide, ...
                 'VerticalAlignment','top', 'HorizontalAlignment','left', 'Interpreter','tex');
        end
        % f_cutoff vertical guide
        if isfinite(pt(a).f_cutoff)
            xline(pt(a).f_cutoff, '--', 'Color', guideCol, 'LineWidth', 1.2, ...
                  'HandleVisibility','off');
        end
    end
    yline(0, '-', 'Color', 0.7*[1 1 1], 'HandleVisibility','off');
    set(gca,'XScale','log','FontSize',fs.axis); box off
    if ~isempty(xLimHz), xlim(xLimHz); end
    xlabel('Flicker frequency (Hz)','FontSize',fs.label);
    ylabel(primLab,'FontSize',fs.label);
    if ismember(lower(groupMode), lower(titleModes))
        title(sprintf('PPC following curve by %s', groupLab), 'FontSize', fs.title);
    end
    legend(h, legLab, 'Location', legendLoc, 'FontSize', fs.legend, 'Box', 'off');

    figName = sprintf('following_descriptors_PPC_%s_%s', groupMode, bootUnit);
    exportFigure(gcf, fullfile(strFigDir, [figName '.png']), exportDPI);
    fprintf('Figure saved: %s\n', fullfile(strFigDir, [figName '.png']));
end

fprintf('\nCaption key (paste into the figure legend):\n  %s\n', keyLine);

fprintf('\nDone. %d PPC descriptor figures written to %s\n', numel(modeList), strFigDir);

%% ####################################################################
%% ##  HELPERS  (verbatim from HKS_following_descriptorsV6.m)
%% ####################################################################

% ---- curve geometry: f_max + f_half on the descending flank --------
% Shared helper: duplicated verbatim. Edit every copy or none (see README).
function d = extractCurveDescriptors(freqList, mu, halfFrac)
    d = struct('f_max',NaN,'f_half',NaN,'peak',NaN,'kMax',NaN,'halfCensored',false);
    f = freqList(:)'; mu = mu(:)';
    idx = find(isfinite(mu));
    if isempty(idx), return; end
    ff = f(idx); mm = mu(idx);
    [peak, kk] = max(mm);
    d.peak = peak;
    d.kMax = idx(kk);                 % index into the FULL grid
    if ~isfinite(peak) || peak <= 0   % no positive following -> descriptors undefined
        return;
    end
    d.f_max = ff(kk);
    thr = halfFrac * peak;            % amplitude: 0.5*peak ; power3dB: peak/sqrt(2)
    for j = kk:(numel(mm)-1)
        a = mm(j); b = mm(j+1);
        if a >= thr && b < thr
            d.f_half = ff(j) + (thr - a)/(b - a) * (ff(j+1) - ff(j));  % linear interp
            return;
        end
    end
    d.halfCensored = true;            % never fell below half within tested range
end

% ---- f_cutoff: significance crossing vs the pre-trial blank --------
% Shared helper: duplicated verbatim. Edit every copy or none (see README).
function [f_cut, sig, pRaw, pAdj, censored] = cutoffCrossing( ...
        freqList, D, kMax, alpha, nLostRun, minCells)
    nFreq = numel(freqList);
    pRaw  = nan(1,nFreq);  posMed = false(1,nFreq);
    for k = 1:nFreq
        x = D(:,k); x = x(isfinite(x));
        if numel(x) >= minCells && any(x ~= 0)
            posMed(k) = median(x) > 0;
            try
                pRaw(k) = signrank(x, 0, 'tail', 'right');   % one-sided: median > 0
            catch
            end
        end
    end
    pAdj = bhfdr(pRaw);
    sig  = (pAdj < alpha) & posMed;

    f_cut = NaN; censored = false;
    if isnan(kMax) || ~any(sig), return; end

    runStart = NaN;
    for k = kMax:(nFreq - nLostRun + 1)
        if all(~sig(k:k+nLostRun-1)), runStart = k; break; end
    end
    if isnan(runStart)
        lastSig  = find(sig, 1, 'last');
        f_cut    = freqList(lastSig);
        censored = true;
    else
        lastSig = find(sig(1:runStart-1), 1, 'last');
        if ~isempty(lastSig) && lastSig >= kMax
            f_cut = freqList(lastSig);
        end
    end
end

% ---- faster bootstrap cutoff: curve crossing vs a fixed floor ------
% Shared helper: duplicated verbatim. Edit every copy or none (see README).
function f_cut = crossingByFloor(freqList, mu, kMax, floorLevel, nLostRun)
    f_cut = NaN;
    if isnan(kMax) || ~isfinite(floorLevel), return; end
    above = mu > floorLevel;
    runStart = NaN;
    for k = kMax:(numel(mu) - nLostRun + 1)
        if all(~above(k:k+nLostRun-1)), runStart = k; break; end
    end
    if isnan(runStart)
        lastA = find(above, 1, 'last');
    else
        lastA = find(above(1:runStart-1), 1, 'last');
    end
    if ~isempty(lastA) && lastA >= kMax, f_cut = freqList(lastA); end
end

% ---- Benjamini-Hochberg FDR (NaN-safe) -----------------------------
% Shared helper: duplicated verbatim. Edit every copy or none (see README).
function padj = bhfdr(p)
    padj = nan(size(p));
    idx  = find(isfinite(p));
    if isempty(idx), return; end
    [ps, ord] = sort(p(idx));
    m   = numel(ps);
    adj = ps(:)' .* m ./ (1:m);
    adj = fliplr(cummin(fliplr(adj)));   % enforce monotonic
    adj = min(adj, 1);
    padj(idx(ord)) = adj;
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
function exportFigure(figHandle, savePath, dpi)
    try
        exportgraphics(figHandle, savePath, 'Resolution', dpi, ...
                       'BackgroundColor', 'white');
    catch
        saveas(figHandle, savePath);
    end
end

% ---- tiny inline ternary -------------------------------------------
% Shared helper: duplicated verbatim. Edit every copy or none (see README).
function s = tern(cond, a, b)
    if cond, s = a; else, s = b; end
end

