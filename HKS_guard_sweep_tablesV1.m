%% HKS_guard_sweep_tablesV1.m
% Spike-guard sensitivity sweep -> supplement tables (NO figures).
% Downstream consumer of HKS_pool_resultantsV3.m -- NO SPIKES TOUCHED.
%
% The skeleton, metric formulas, descriptor logic, cutoff test and hierarchical
% bootstrap match HKS_thesis_descriptors_PPCV1.m: the same Mcurve/Dprim split
% (raw PPC for the curve, paired difference for the cutoff) and the same null
% gate. That is what makes any number here directly comparable to the
% main-analysis run at minSpkPerBin = 10. Nothing about the estimator changes;
% only the guard moves.
%
% ---------------------------------------------------------------------
% INPUT   <DRIVE>\Kersten\Pooled\resonance_resultants_pool.mat
% OUTPUT  <DRIVE>\Kersten\Tables\guard_sweep\guardsweep_<metric>_<date>.*
% ---------------------------------------------------------------------
%   T_desc      one row per (mode, group, guard): point estimate + boot CI for
%               all three descriptors, plus n's and censor/gate flags. Raw
%               ungated values are kept in the *_ungated columns so the gate is
%               inspectable, not silent.
%   T_contrast  one row per (mode, guard, descriptor, groupA, groupB): paired
%               within-iteration median difference + CI + excludes-zero
%   T_surv      one row per (mode, group, guard, frequency): paired-valid cell
%               count -- the guard-survival / selection-bias audit
%   T_stab      one row per (mode, group, descriptor): spread of the point
%               estimate across the whole guard set + deviation from the
%               reference guard. This is the table the one-line insensitivity
%               statement in the methods is actually claiming.
%   Written as CSV, one XLSX workbook (one sheet per table), and a .txt of the
%   formatted groups x guards matrices for direct paste into the supplement.
%
% ---------------------------------------------------------------------
% CURRENT BEHAVIOUR
% ---------------------------------------------------------------------
%   * ONE global guard per run, applied identically to every group, so the
%     between-group contrasts stay interpretable. No per-group guards.
%   * The guard is set by hand and the sweep LICENSES that choice. Deriving a
%     threshold from the data would reintroduce the circularity the sweep exists
%     to remove.
%   * Inclusion, bootstrap unit, estimator and seed are frozen across the sweep.
%     One knob moves; anything else changing makes the comparison meaningless.
%   * rng is re-seeded per (mode, guard), so every cell of the sweep draws the
%     same resample stream and differences between cells are the guard, not the
%     RNG.
%   * NULL GATE. f_max and f_half are pure curve geometry, so max() of a noise
%     curve always returns SOME frequency and a group with no real following
%     (CA1) would get finite-looking descriptors with absurd CIs. When no
%     frequency survives the paired test, all three descriptors go NaN and the
%     row is marked n.d.
%   * No figures. The sweep's job is the flatness claim, not a display.
%
% ---------------------------------------------------------------------
% RUN ORDER    HKS_pool_resultantsV3 -> THIS  (supplement to
%              HKS_thesis_descriptors_PPCV1, which it must match)
% DEPENDS ON   Statistics Toolbox (signrank, prctile)
% ---------------------------------------------------------------------

clear all; close all  %#ok<CLALL>

%% ####################################################################
%% ##  SETUP
%% ####################################################################

%% ---- [S1] Pooled input file + output dir -------------------------
driveLetter = 'W:';                     % drive letter for the data root
poolFile    = fullfile(driveLetter, 'Kersten', 'Pooled', ...
                       'resonance_resultants_pool.mat');
outDir      = fullfile(driveLetter, 'Kersten', 'Tables', 'guard_sweep');

writeCSV  = true;    % one CSV per table
writeXLSX = true;    % one workbook, one sheet per table (nicer for a supplement)
writeTXT  = true;    % the formatted groups x guards matrices, as printed

%% ---- [S2] Inclusion mask (must match the descriptor script) ------
% Guard-independent: the sweep must move ONE knob, so inclusion is frozen here.
useKS         = false;
useBC         = true;
useArtefact   = false;
useResponsive = true;

%% ---- [S3] The sweep itself ---------------------------------------
minSpkList = [5 10 20 50];              % guard values to sweep (spikes/bin)
refGuard   = 10;                        % the value reported in the thesis body
modeList   = {'area', 'v1layer', 'v1celltype'};
% add 'v1layertype' below if the supplement wants the thin layer x type cross;
% expect small-n warnings and wide CIs there.
% modeList = {'area','v1layer','v1celltype','v1layertype'};

%% ---- [S4] Group definitions (must match the descriptor script) ---
layerDefs = { 'L2/3', {"L2/3"} ; ...
              'L4',   {"L4"}   ; ...
              'L5',   {"L5"}   ; ...
              'L6',   {"L6a","L6b"} };
typeDefs   = { 'NS', {"NS"} ; ...
               'BS', {"BS"} };
NS_CUTOFF_US = 400;          % peak-trough; < = NS, >= = BS (fallback derivation only)

%% ---- [S5] Metric + curve parameters ------------------------------
metricPrimary  = 'PPC';      % descriptors read from here (locked primary)
absUnit        = 'cycle';    % |F1| per-cycle normalisation; live only when
                             % metricPrimary = 'F1' (not swept)
halfMode       = 'amplitude';   % 'amplitude' -> 0.5*peak ; 'power3dB' -> peak/sqrt(2)
switch lower(halfMode)
    case 'amplitude', halfFrac = 0.5;
    case 'power3db',  halfFrac = 1/sqrt(2);
    otherwise, error('halfMode must be amplitude | power3dB');
end

%% ---- [S6] f_cutoff test parameters -------------------------------
cutoffAlpha     = 0.05;      % BH-FDR level across the frequency grid
nLostRun        = 2;         % consecutive n.s. bins that define "staying lost"
minCellsForTest = 5;         % min finite paired cells to attempt signrank at a freq

%% ---- [S7] Bootstrap parameters -----------------------------------
doBootstrap = true;
bootUnit    = 'animal';      % 'animal' (n=5, honest) | 'recording'
nBoot       = 1000;          % numel(modeList) x numel(minSpkList) bootstraps at
                             % this many iters is the runtime driver
bootCI      = [2.5 97.5];
cutoffBootMethod = 'signrank';   % 'signrank' = faithful to the main run (SLOW).
                                 % 'curveFloor' is ~10-30x faster but is a
                                 % DIFFERENT estimator -- only use it if the
                                 % main run also used it, else the supplement
                                 % is not comparable to the body.
rngSeed = 7;                 % re-seeded per (mode,guard) so cells are comparable

%% ---- [S8] Null gate + reporting ----------------------------------
gateOnSignificance = true;   % no significant freq -> all descriptors NaN (n.d.)
prettyUses = 'point';        % 'point' (point estimate + boot CI) | 'bootmed'
ndMark     = 'n.d.';         % printed for gated groups

fprintf('\n=== Guard sweep : %s primary, guards [%s] ===\n', ...
        metricPrimary, num2str(minSpkList));

%% ---- Load pooled file --------------------------------------------
assert(exist(poolFile,'file')==2, 'Pooled file not found: %s', poolFile);
P = load(poolFile);
C = P.clusterTable;
freqList = P.freqList(:)';
nFreq    = numel(freqList);
fprintf('Loaded %d clusters x %d freqs from %s\n', height(C), nFreq, poolFile);

%% ---- Inclusion mask ----------------------------------------------
incl = true(height(C),1);
if useKS,         incl = incl & C.pass_ksLabel;     end
if useBC,         incl = incl & C.pass_bcLabel;     end
if useArtefact,   incl = incl & C.pass_notArtefact; end
if useResponsive, incl = incl & C.pass_responsive;  end
fprintf('Inclusion (KS=%d BC=%d art=%d resp=%d): %d / %d clusters.\n', ...
        useKS,useBC,useArtefact,useResponsive, sum(incl), height(C));

region  = strtrim(string(C.region));
mouseID = C.mouseID;
recID   = strtrim(string(C.recID));

%% ---- Layer / cell-type labels (same parse as the descriptor script) ---
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
    fprintf('v1Layer: derived in place from regionAllen.\n');
end

cellType = repmat("", height(C), 1);
needCellType = any(ismember(lower(modeList), {'v1celltype','v1layertype'}));
if needCellType
    if any(strcmp(C.Properties.VariableNames, 'cellType'))
        cellType = strtrim(string(C.cellType));
        fprintf('cellType: using pooled column.\n');
    elseif any(strcmp(C.Properties.VariableNames, 'wfDuration_us'))
        wf = double(C.wfDuration_us);
        cellType(wf <  NS_CUTOFF_US) = "NS";
        cellType(wf >= NS_CUTOFF_US) = "BS";
        fprintf('cellType: derived in place from wfDuration_us at %g us.\n', NS_CUTOFF_US);
    else
        error(['Pooled table has neither cellType nor wfDuration_us. Rerun: ' ...
               'HKS_bombcell_batchV3 -> HKS_cluster_QC_unifiedV7 -> ' ...
               'HKS_pool_resultantsV3.']);
    end
end

%% ####################################################################
%% ##  [1] COMPUTE : guard-INDEPENDENT metrics (computed once)
%% ####################################################################
% Metric formulas verbatim from HKS_fourierV4.m D.1-D.3. Only the pairOK mask
% below depends on the guard, so everything up to the mask is computed a single
% time and reused across the sweep.
Zc = P.Zc; Nc = P.Nc; Tc = P.Tc;
Zn = P.Zn; Nn = P.Nn; Tn = P.Tn;

VS_r  = abs(Zc);
PPC_r = (Nc .* VS_r.^2 - 1) ./ max(Nc - 1, 1);           % Vinck 2010
F1_r  = 2 .* Nc .* VS_r ./ max(Tc, eps);                 % spikes/s
VS_n  = abs(Zn);
PPC_n = (Nn .* VS_n.^2 - 1) ./ max(Nn - 1, 1);
F1_n  = 2 .* Nn .* VS_n ./ max(Tn, eps);

if strcmpi(absUnit,'cycle')
    F1_r = F1_r ./ freqList;   F1_n = F1_n ./ freqList;
end

D_PPC_raw = PPC_r - PPC_n;      % unguarded; masked per guard inside the loop
D_F1_raw  = F1_r  - F1_n;

switch upper(metricPrimary)
    case 'PPC', Mraw = PPC_r;  Draw = D_PPC_raw;
    case 'F1',  Mraw = F1_r;   Draw = D_F1_raw;
    otherwise, error('metricPrimary must be PPC | F1');
end

%% ---- bootstrap unit list (guard-independent, so CIs stay comparable)
switch lower(bootUnit)
    case 'animal',    unitKey = mouseID;  unitList = unique(mouseID(incl));
    case 'recording', unitKey = recID;    unitList = unique(recID(incl));
    otherwise, error('bootUnit must be animal | recording');
end
nUnit = numel(unitList);
fprintf('Bootstrap unit: %s, n=%d units, %d iters, method=%s\n', ...
        bootUnit, nUnit, nBoot, cutoffBootMethod);
if nUnit < 4
    warning('Only %d bootstrap units -- CIs will be very coarse.', nUnit);
end

%% ####################################################################
%% ##  [2] SWEEP : groupMode x minSpkPerBin
%% ####################################################################
% Preallocated with the right COLUMN count: `X(end+1,:) = {...}` on a 0x0 cell
% errors (LHS is 1x0), so the row width has to exist up front.
descRows = cell(0,27);  contRows = cell(0,10);  survRows = cell(0,5);
prettyTxt = strings(0,1);                        % formatted matrices for .txt
nGuard = numel(minSpkList);
tSweep = tic;

for iM = 1:numel(modeList)
    groupMode = modeList{iM};

    % ---- build the group set for this mode (masks, not handles) ----
    [grpLab, grpMask, groupLab] = buildGroupSet(groupMode, region, v1Layer, ...
                                        cellType, layerDefs, typeDefs);
    nGrp = numel(grpLab);
    grpRows = cell(nGrp,1);
    for g = 1:nGrp, grpRows{g} = find(incl & grpMask{g}); end

    fprintf('\n###############################################\n');
    fprintf('## MODE: %s  (%s, %d groups)\n', groupMode, groupLab, nGrp);
    fprintf('###############################################\n');
    for g = 1:nGrp
        fprintf('   %-10s  n=%4d cells, %d animals\n', grpLab{g}, ...
                numel(grpRows{g}), numel(unique(mouseID(grpRows{g}))));
    end
    if any(cellfun(@isempty, grpRows))
        warning('Mode %s: group(s) with zero included cells -- rows will be NaN.', groupMode);
    end

    % storage for the pretty matrices: {descriptor}(group, guard) strings
    prettyCell = repmat({repmat("", nGrp, nGuard)}, 1, 3);
    descNames  = {'f_max','f_half','f_cutoff'};

    for iG = 1:nGuard
        minSpkPerBin = minSpkList(iG);
        tGuard = tic;
        rng(rngSeed, 'twister');     % identical resample stream per cell of the sweep

        % ---- [2.1] apply the guard ------------------------------------
        realOK = (Nc >= minSpkPerBin);
        pairOK = realOK & (Nn >= minSpkPerBin);
        Mcurve = Mraw;  Mcurve(~realOK) = NaN;   % f_max / f_half / kMax
        Dprim  = Draw;  Dprim(~pairOK) = NaN;    % paired test -> f_cutoff

        fprintf('\n-- %s | guard = %d spikes/bin | paired-valid %d/%d entries --\n', ...
                groupMode, minSpkPerBin, sum(pairOK(:)), numel(pairOK));

        % ---- [2.2] point-estimate descriptors + survival --------------
        ptF   = nan(nGrp,3);   ptRawF = nan(nGrp,3);
        sigAny = false(nGrp,1);
        halfCens = false(nGrp,1);  cutCens = false(nGrp,1);
        peakV = nan(nGrp,1);
        nValid = nan(nGrp, nFreq);

        for g = 1:nGrp
            rows = grpRows{g};
            if isempty(rows), continue; end
            mu = mean(Mcurve(rows,:), 1, 'omitnan');
            nValid(g,:) = sum(isfinite(Dprim(rows,:)), 1);

            cd = extractCurveDescriptors(freqList, mu, halfFrac);
            [f_cut, sigVec, ~, ~, cCens] = cutoffCrossing( ...
                freqList, Dprim(rows,:), cd.kMax, cutoffAlpha, nLostRun, minCellsForTest);

            ptRawF(g,:) = [cd.f_max, cd.f_half, f_cut];
            sigAny(g)   = any(sigVec);
            halfCens(g) = cd.halfCensored;
            cutCens(g)  = cCens;
            peakV(g)    = cd.peak;

            % null gate: no frequency survives the paired test -> n.d.
            if gateOnSignificance && ~sigAny(g)
                ptF(g,:) = NaN;
            else
                ptF(g,:) = ptRawF(g,:);
            end

            fprintf('   %-10s n=%4d | f_max=%6.1f f_half=%6.1f%s f_cutoff=%6.1f%s %s\n', ...
                grpLab{g}, numel(rows), ptF(g,1), ptF(g,2), tern(halfCens(g),'*',' '), ...
                ptF(g,3), tern(cutCens(g),'*',' '), tern(~sigAny(g), '<- gated (n.d.)',''));
        end

        % ---- [2.3] hierarchical bootstrap -----------------------------
        bLo = nan(nGrp,3); bHi = nan(nGrp,3); bMed = nan(nGrp,3);
        bootRaw = cell(1,3);
        if doBootstrap
            unitRowsPerGrp = cell(nGrp, nUnit);
            for g = 1:nGrp
                for u = 1:nUnit
                    unitRowsPerGrp{g,u} = grpRows{g}(unitKey(grpRows{g}) == unitList(u));
                end
            end
            floorLevel = nan(nGrp,1);
            for g = 1:nGrp
                if isfinite(ptRawF(g,3))
                    [~,kc] = min(abs(freqList - ptRawF(g,3)));
                    % Mcurve, not Dprim: crossingByFloor compares this level
                    % against mu, and mu is the raw-PPC curve.
                    floorLevel(g) = mean(Mcurve(grpRows{g},kc), 1, 'omitnan');
                end
            end

            bM = nan(nBoot,nGrp); bH = nan(nBoot,nGrp); bC = nan(nBoot,nGrp);
            for b = 1:nBoot
                draw = unitList(randi(nUnit, nUnit, 1));
                for g = 1:nGrp
                    sel = cell(numel(draw),1);
                    for d = 1:numel(draw)
                        uu = find(unitList == draw(d), 1);
                        sel{d} = unitRowsPerGrp{g,uu};
                    end
                    rows = vertcat(sel{:});
                    if isempty(rows), continue; end
                    mu = mean(Mcurve(rows,:), 1, 'omitnan');
                    cd = extractCurveDescriptors(freqList, mu, halfFrac);
                    bM(b,g) = cd.f_max;
                    bH(b,g) = cd.f_half;
                    switch lower(cutoffBootMethod)
                        case 'signrank'
                            bC(b,g) = cutoffCrossing(freqList, Dprim(rows,:), cd.kMax, ...
                                                     cutoffAlpha, nLostRun, minCellsForTest);
                        case 'curvefloor'
                            bC(b,g) = crossingByFloor(freqList, mu, cd.kMax, ...
                                                      floorLevel(g), nLostRun);
                        otherwise, error('cutoffBootMethod must be signrank | curveFloor');
                    end
                end
            end
            % Gate the CIs on the same rule as ptF, or a gated group still
            % returns an interval and every contrast against it compares with
            % an undefined quantity.
            if gateOnSignificance
                for g = 1:nGrp
                    if ~sigAny(g)
                        bM(:,g) = NaN;  bH(:,g) = NaN;  bC(:,g) = NaN;
                    end
                end
            end

            bootRaw = {bM, bH, bC};
            for k = 1:3
                M = bootRaw{k};
                bLo(:,k)  = prctile(M, bootCI(1), 1)';
                bHi(:,k)  = prctile(M, bootCI(2), 1)';
                bMed(:,k) = median(M, 1, 'omitnan')';
            end

            % ---- [2.4] paired between-group contrasts ------------------
            if nGrp >= 2
                pairs = nchoosek(1:nGrp, 2);
                for k = 1:3
                    M = bootRaw{k};
                    for pI = 1:size(pairs,1)
                        i = pairs(pI,1); j = pairs(pI,2);
                        dd = M(:,i) - M(:,j);
                        lo = prctile(dd, bootCI(1));
                        hi = prctile(dd, bootCI(2));
                        md = median(dd,'omitnan');
                        gated = gateOnSignificance && (~sigAny(i) || ~sigAny(j));
                        contRows(end+1,:) = { string(groupMode), minSpkPerBin, ...
                                string(descNames{k}), string(grpLab{i}), string(grpLab{j}), ...
                                md, lo, hi, ~(lo<=0 && hi>=0), gated };  %#ok<SAGROW>
                    end
                end
            end
        end

        % ---- [2.5] rows for the flat tables ---------------------------
        for g = 1:nGrp
            rows = grpRows{g};
            useV = ptF(g,:);
            if strcmpi(prettyUses,'bootmed')
                pv = bMed(g,:);
                if gateOnSignificance && ~sigAny(g), pv = nan(1,3); end
            else
                pv = useV;
            end
            for k = 1:3
                cens = false;
                if k==2, cens = halfCens(g); elseif k==3, cens = cutCens(g); end
                prettyCell{k}(g,iG) = fmtCI(pv(k), bLo(g,k), bHi(g,k), cens, ...
                                            sigAny(g) || ~gateOnSignificance, ndMark);
            end

            descRows(end+1,:) = { string(groupMode), string(grpLab{g}), minSpkPerBin, ...
                numel(rows), numel(unique(mouseID(rows))), ...
                median(nValid(g,:),'omitnan'), min(nValid(g,:)), max(nValid(g,:)), ...
                sigAny(g), peakV(g), ...
                ptF(g,1), bLo(g,1), bHi(g,1), bMed(g,1), ...
                ptF(g,2), bLo(g,2), bHi(g,2), bMed(g,2), halfCens(g), ...
                ptF(g,3), bLo(g,3), bHi(g,3), bMed(g,3), cutCens(g), ...
                ptRawF(g,1), ptRawF(g,2), ptRawF(g,3) };  %#ok<SAGROW>

            for k = 1:nFreq
                survRows(end+1,:) = { string(groupMode), string(grpLab{g}), ...
                    minSpkPerBin, freqList(k), nValid(g,k) };  %#ok<SAGROW>
            end
        end

        fprintf('   [%s | guard %d done in %.1f s]\n', groupMode, minSpkPerBin, toc(tGuard));
    end

    % ---- [2.6] pretty groups x guards matrices for this mode ---------
    for k = 1:3
        blk = prettyBlock(sprintf('%s  --  %s  (%s)', upper(groupMode), ...
                    descNames{k}, groupLab), grpLab, minSpkList, prettyCell{k}, refGuard);
        fprintf('%s', blk);
        prettyTxt(end+1,1) = string(blk); %#ok<SAGROW>
    end
end
fprintf('\nSweep finished in %.1f min.\n', toc(tSweep)/60);

%% ####################################################################
%% ##  [3] ASSEMBLE tables
%% ####################################################################
T_desc = cell2table(descRows, 'VariableNames', { ...
    'mode','group','guard','nCell','nAnimal', ...
    'nValid_med','nValid_min','nValid_max','anySigFreq','peak', ...
    'f_max','f_max_lo','f_max_hi','f_max_bmed', ...
    'f_half','f_half_lo','f_half_hi','f_half_bmed','f_half_censored', ...
    'f_cutoff','f_cut_lo','f_cut_hi','f_cut_bmed','f_cut_censored', ...
    'f_max_ungated','f_half_ungated','f_cut_ungated'});

if ~isempty(contRows)
    T_contrast = cell2table(contRows, 'VariableNames', { ...
        'mode','guard','descriptor','groupA','groupB', ...
        'medianDiff','lo','hi','excludesZero','gated'});
else
    T_contrast = table();
end

T_surv = cell2table(survRows, 'VariableNames', { ...
    'mode','group','guard','freq_Hz','nPairedValid'});

%% ---- [3b] stability summary (the claim the methods line makes) -----
% For each (mode, group, descriptor): min/max/range of the POINT estimate over
% the swept guards, and the largest deviation from the reference guard. A
% descriptor is "insensitive" when range is within a couple of grid steps and
% the sign of every between-group contrast is unchanged (see T_contrast).
stabRows = cell(0,10);
descNames = {'f_max','f_half','f_cutoff'};
keys = unique(T_desc(:,{'mode','group'}), 'stable');   % tables: row-wise by default
for r = 1:height(keys)
    sub = T_desc(T_desc.mode==keys.mode(r) & T_desc.group==keys.group(r), :);
    for k = 1:3
        v  = sub.(descNames{k});
        gv = sub.guard;
        vRef = v(gv==refGuard);
        if isempty(vRef), vRef = NaN; else, vRef = vRef(1); end
        vf = v(isfinite(v));
        if isempty(vf)
            mn = NaN; mx = NaN; rg = NaN; md = NaN;
        else
            mn = min(vf); mx = max(vf); rg = mx - mn;
            md = max(abs(v - vRef), [], 'omitnan');
        end
        stabRows(end+1,:) = { keys.mode(r), keys.group(r), string(descNames{k}), ...
            vRef, mn, mx, rg, md, sum(isfinite(v)), numel(v) };  %#ok<SAGROW>
    end
end
T_stab = cell2table(stabRows, 'VariableNames', { ...
    'mode','group','descriptor','value_at_ref','min_over_guards','max_over_guards', ...
    'range_over_guards','maxAbsDev_from_ref','nGuards_defined','nGuards_tested'});

fprintf('\n--- Stability across guards [%s] (point estimates, ref = %d) ---\n', ...
        num2str(minSpkList), refGuard);
disp(T_stab);

% contrast-sign stability: does any pairwise verdict flip across the sweep?
if ~isempty(contRows)
    ck = unique(T_contrast(:,{'mode','descriptor','groupA','groupB'}), 'stable');
    nFlip = 0;
    fprintf('\n--- Contrast verdict stability (excludes-zero across guards) ---\n');
    for r = 1:height(ck)
        m = T_contrast.mode==ck.mode(r) & T_contrast.descriptor==ck.descriptor(r) & ...
            T_contrast.groupA==ck.groupA(r) & T_contrast.groupB==ck.groupB(r);
        ez = T_contrast.excludesZero(m);
        if numel(unique(ez)) > 1
            nFlip = nFlip + 1;
            fprintf('  FLIPS: %-12s %-9s %-8s vs %-8s : [%s]\n', ck.mode(r), ...
                ck.descriptor(r), ck.groupA(r), ck.groupB(r), num2str(ez(:)'));
        end
    end
    if nFlip == 0
        fprintf('  No verdict flips: every pairwise contrast keeps its excludes-zero\n');
        fprintf('  status across all %d guards. This is the strong form of the claim.\n', nGuard);
    else
        fprintf('  %d contrast(s) flip across the sweep -- report these explicitly.\n', nFlip);
    end
end

%% ####################################################################
%% ##  [4] EXPORT
%% ####################################################################
if ~exist(outDir,'dir'), mkdir(outDir); end
stamp = char(datetime('today','Format','yyyy-MM-dd'));
tag   = sprintf('guardsweep_%s_%s', metricPrimary, stamp);

if writeCSV
    writetable(T_desc,     fullfile(outDir, [tag '_descriptors.csv']));
    writetable(T_stab,     fullfile(outDir, [tag '_stability.csv']));
    writetable(T_surv,     fullfile(outDir, [tag '_survival.csv']));
    if ~isempty(contRows)
        writetable(T_contrast, fullfile(outDir, [tag '_contrasts.csv']));
    end
    fprintf('\nCSVs written to %s\n', outDir);
end

if writeXLSX
    xlsFile = fullfile(outDir, [tag '.xlsx']);
    if exist(xlsFile,'file'), delete(xlsFile); end
    writetable(T_desc, xlsFile, 'Sheet', 'descriptors');
    writetable(T_stab, xlsFile, 'Sheet', 'stability');
    if ~isempty(contRows)
        writetable(T_contrast, xlsFile, 'Sheet', 'contrasts');
    end
    writetable(T_surv, xlsFile, 'Sheet', 'survival');
    fprintf('Workbook written: %s\n', xlsFile);
end

if writeTXT
    txtFile = fullfile(outDir, [tag '_tables.txt']);
    fid = fopen(txtFile, 'w');
    fprintf(fid, 'Spike-guard sweep, %s primary, guards [%s], ref = %d\n', ...
            metricPrimary, num2str(minSpkList), refGuard);
    fprintf(fid, 'bootUnit=%s (n=%d), nBoot=%d, cutoffBootMethod=%s, CI=[%g %g]\n', ...
            bootUnit, nUnit, nBoot, cutoffBootMethod, bootCI(1), bootCI(2));
    fprintf(fid, 'Inclusion: KS=%d BC=%d artefact=%d responsive=%d\n', ...
            useKS, useBC, useArtefact, useResponsive);
    fprintf(fid, 'Cells are  point [lo-hi] ; * = right-censored ; %s = gated (no significant frequency)\n\n', ndMark);
    for i = 1:numel(prettyTxt)
        fprintf(fid, '%s', prettyTxt(i));
    end
    fclose(fid);
    fprintf('Formatted tables written: %s\n', txtFile);
end

fprintf(['\nDone. T_desc / T_contrast / T_surv / T_stab are in the workspace.\n' ...
         'For the methods line, read range_over_guards in T_stab and the\n' ...
         'verdict-stability block above.\n']);

%% ####################################################################
%% ##  HELPERS  (local functions; keep at end of file)
%% ####################################################################

% ---- group set as masks (mode -> labels + logical masks) -----------
function [lab, mask, groupLab] = buildGroupSet(groupMode, region, v1Layer, ...
                                        cellType, layerDefs, typeDefs)
    switch lower(groupMode)
        case 'area'
            % display label 'dLGN' maps to the stored region string 'LGN'
            areaMap = { 'dLGN', 'LGN' ; 'V1', 'V1' ; 'CA1', 'CA1' };
            lab  = areaMap(:,1);
            mask = cell(size(lab));
            for a = 1:numel(lab), mask{a} = strcmp(region, areaMap{a,2}); end
            groupLab = 'area';
        case 'v1layer'
            lab  = layerDefs(:,1);
            mask = cell(size(lab));
            for a = 1:numel(lab)
                mask{a} = (region=="V1") & ismember(v1Layer, string(layerDefs{a,2}));
            end
            groupLab = 'V1 layer';
        case 'v1celltype'
            lab  = typeDefs(:,1);
            mask = cell(size(lab));
            for a = 1:numel(lab)
                mask{a} = (region=="V1") & ismember(cellType, string(typeDefs{a,2}));
            end
            groupLab = 'V1 cell type';
        case 'v1layertype'
            nL = size(layerDefs,1); nT = size(typeDefs,1);
            lab = cell(nL*nT,1); mask = cell(nL*nT,1); rr = 0;
            for iL = 1:nL
                for iT = 1:nT
                    rr = rr + 1;
                    lab{rr}  = sprintf('%s-%s', layerDefs{iL,1}, typeDefs{iT,1});
                    mask{rr} = (region=="V1") ...
                             & ismember(v1Layer,  string(layerDefs{iL,2})) ...
                             & ismember(cellType, string(typeDefs{iT,2}));
                end
            end
            groupLab = 'V1 layer x cell type';
        otherwise
            error('groupMode must be area | v1layer | v1celltype | v1layertype');
    end
end

% ---- format one table cell:  point [lo-hi] -------------------------
function s = fmtCI(v, lo, hi, censored, ok, ndMark)
    if ~ok || ~isfinite(v)
        s = string(ndMark);
        return;
    end
    if isfinite(lo) && isfinite(hi)
        s = string(sprintf('%.1f [%.1f-%.1f]%s', v, lo, hi, tern(censored,'*','')));
    else
        s = string(sprintf('%.1f%s', v, tern(censored,'*','')));
    end
end

% ---- build one printable groups x guards block ---------------------
function blk = prettyBlock(titleStr, grpLab, guards, cells, refGuard)
    w = 22;
    L = sprintf('\n%s\n', titleStr);
    L = [L sprintf('%s\n', repmat('-', 1, 12 + w*numel(guards)))];
    L = [L sprintf('%-12s', 'group')];
    for iG = 1:numel(guards)
        h = sprintf('guard=%d%s', guards(iG), tern(guards(iG)==refGuard,' (ref)',''));
        L = [L sprintf('%-*s', w, h)]; %#ok<AGROW>
    end
    L = [L sprintf('\n')];
    for g = 1:numel(grpLab)
        L = [L sprintf('%-12s', grpLab{g})]; %#ok<AGROW>
        for iG = 1:numel(guards)
            L = [L sprintf('%-*s', w, cells(g,iG))]; %#ok<AGROW>
        end
        L = [L sprintf('\n')]; %#ok<AGROW>
    end
    blk = L;
end

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
    d.kMax = idx(kk);
    if ~isfinite(peak) || peak <= 0
        return;
    end
    d.f_max = ff(kk);
    thr = halfFrac * peak;
    for j = kk:(numel(mm)-1)
        a = mm(j); b = mm(j+1);
        if a >= thr && b < thr
            d.f_half = ff(j) + (thr - a)/(b - a) * (ff(j+1) - ff(j));
            return;
        end
    end
    d.halfCensored = true;
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
                pRaw(k) = signrank(x, 0, 'tail', 'right');
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
    adj = fliplr(cummin(fliplr(adj)));
    adj = min(adj, 1);
    padj(idx(ord)) = adj;
end

% ---- tiny inline ternary -------------------------------------------
% Shared helper: duplicated verbatim. Edit every copy or none (see README).
function s = tern(cond, a, b)
    if cond, s = a; else, s = b; end
end
