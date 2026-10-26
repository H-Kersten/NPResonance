# Analysis pipeline — README

Single-neuron temporal following to visual flicker in awake head-fixed mice.
Neuropixels 1.0 (SpikeGLX) across **dLGN, V1, CA1**; 48-frequency flicker protocol
(**2–243 Hz**) plus pulse and white-noise blocks. The question is whether the
frequency-banded resonance seen at the population level (EEG / ssVEP) is carried by
individual neurons, and how temporal fidelity degrades along dLGN → V1 → CA1.

Analysis is **MATLAB**. Upstream preprocessing (CatGT → Kilosort4 → Phy2, plus Acquipix
for stream sync and stimulus reconstruction) is already done per recording.

**Cohort:** 5 mice (M5, M6, M7, M8, M10), 356 QC-passing units (dLGN 166, V1 140, CA1 50).
M9 was dropped for probe drift.

---

## Run order

| # | Stage | Script | Output |
|---|-------|--------|--------|
| 1 | Bombcell QC | `HKS_bombcell_batchV3.m` | per-rec `kilosort4/cluster_bc_unitType.tsv` + `bombcell_pooled_QC.mat` |
| 2 | Unified QC | `HKS_cluster_QC_unifiedV7.m` | per-rec `<recName>_cluster_QC.mat` |
| 3 | Stage-1 reduction | `HKS_pool_resultantsV3.m` | `resonance_resultants_pool.mat` |
| 4a | Metrics figures | `HKS_thesis_metrics_pooledV3.m` (+ `_supplementV1`) | Figure 3 panels |
| 4b | Descriptors | `HKS_thesis_descriptors_PPCV2.m` | Figure 4 panels |
| 4c | Exposure figure | `HKS_spike_count_figureV2.m` | Figure 2 panels A–D |
| 4c | Rasters | `HKS_unit_rastersV3.m` | Figure 2 panels E–G |
| 4d | Tables | `HKS_qc_supplement_tableV1.m`, `HKS_sample_tableV1.m`, `HKS_guard_sweep_tablesV1.m` | Results / supplement tables |

**After stage 3 the analysis never touches spikes again.** Everything downstream reduces
the pooled resultant file. Stages 1–3 are the only expensive reruns; stage 4 is cheap.
The one exception is `HKS_unit_rastersV3`, which reads per-recording QC files directly
because it needs individual spike times.

⚠ **Check your file versions first.** Older versions of several scripts are still on the
drive, and the three stage-4 scripts above (`pooledV3`, `PPCV2`, `rastersV3`) were **not
on the Drive at all** at the time of this cleanup — they came off the office PC. Make sure
the Drive carries them before relying on it as the handover copy. Note also that
`thesis_metrics_pooled` jumps V1 → V3 with no V2 in either location. The run order above names the current ones. Running stage 3 with
`pool_resultantsV2` uses the old analysis window and will not reproduce the reported
numbers; running stage 4b with `following_descriptorsV6` will not either (see below).

Paths: data root `W:\Gulbinaite\Resonance_neuropixels\Data_collection\Experiment_1\M<n>\Ephys`;
packages `W:\Kersten\Packages_and_Resources`; figures `W:\Kersten\Figures\...`.

---

## Conventions

- **Naming.** `RMS_*` = Rasa's originals, read-only reference. `HKS_*` = working scripts,
  versioned `…V1`, `…V2`. Headers describe **current behaviour**, not deltas from a
  previous version.
- **Layout.** `SETUP → COMPUTE → PLOTS`, with `[S#]` / `[#]` / `[P-#]` block numbering —
  all maths finishes before any `figure()` call. `fullfile` everywhere, no `cd`, no
  trailing slashes. `-v7.3` MAT saves. Comments are ASCII-only (`Delta`, `us`, `->`).
- **Recording selection.** Batch scripts share a `recSpec = {mouseID, 'all' | {dates}}`
  format. Folders are found by globbing `catgt*` and parsing the date. There is one
  recording per day, so **the date is the unique key** — the `R01` token is not a run
  index. Each `recID` must resolve to **exactly one** `catgt*` folder; a duplicate on disk
  is a hard error rather than a silent choice.
- **Keys.** `globalID = "<recName>__clu<id>"` (cross-recording unique).
  `recID = "<date>R0N"`. Cluster handle is `sSpikes.cids`, 0-indexed, matches Phy.
- **Store continuous, threshold later.** QC metrics, spike guards and inclusion criteria
  are stored as continuous values; thresholds are applied at load time. Adding an audit
  column is always preferred to hard-coding a cut.

---

## Stage 1 — Bombcell QC

### `HKS_bombcell_batchV3.m`

**In:** raw `*.ap.bin` + `.meta`, `kilosort4/` folder, `recSpec`.
**Out:** per-recording `<ksdir>/cluster_bc_unitType.tsv` (columns `cluster_id`,
`bc_unitType`) — this is the handoff to stage 2 — plus pooled
`W:\Kersten\QC\bombcell_pooled\bombcell_pooled_QC.mat`.

Runs Bombcell per recording and pools every unit into one table keyed by `globalID`.
Uses string unit-type labels with the **non-somatic split ON**, so `NOISE` and
`NON-SOMA` stay distinguishable; `excludeLabels = {'NOISE','NON-SOMA'}` keeps GOOD plus
somatic MUA.

**Custom RPV promotion.** Each recording is classified twice — once at the default
refractory-period-violation threshold, once relaxed — and somatic **MUA is promoted to
GOOD** where the unit was MUA *solely* on RPV. Applied before the TSV is written.
`NON-SOMA` and `NOISE` are never promoted. Original label kept as `bcLabelRaw`, the
promotion flagged in `rpvPromoted`. At `rpvPromoteCutoff = 0.15` this promoted 124 of
3,005 somatic MUA on top of 696 native GOOD; the promoted units are statistically
indistinguishable from native GOOD on SNR and presence ratio.

⚠ **`rpvPromoteCutoff = 0.15` is provisional and unvalidated.** It sets `bc_unitType`,
hence `pass_bcLabel`, hence the analysed cohort. The value was chosen by hand: the
25th-percentile route gave ~0.08, which left too few units to work with, and 0.1 is
regarded as strict in the lab. **No sensitivity analysis backs it**, and it was chosen
partly on yield — a researcher degree of freedom worth stating plainly rather than
discovering later. Treat it as an open parameter, and ideally sweep it the way
`minSpkPerBin` is swept: `passMatBC` and `bcLabelRaw` are already saved in
`bombcell_pooled_QC.mat`, so a sweep can run off the saved pool without re-running the
heavy Bombcell loop.

⚠ **The percentile gate never reaches an analysis.** `inclPctile = 25`, `incl_pctl` and
`keep` are computed and stored, but the TSV carries `bc_unitType` only, and downstream QC
sets `pass_bcLabel` from that label. The percentile columns are a descriptive
distribution, not an exclusion step — **don't draw them as a stage of the inclusion
funnel** in any figure or supplement table.

**Two accountings, kept separate on purpose.** `metricPassPctl` (percentile, descriptive)
and `metricPassBC` (Bombcell thresholds, gating). In both, `nOnlyFail` is a metric's
*marginal* cost and is the honest number to quote; the plain `nFail` counts overlap
heavily and sum to nothing meaningful. The label-reconstruction audit quantifies how much
of the GOOD call the 9 `metricSpec` metrics account for — Bombcell also applies
waveform-shape and spikes-missing criteria that are not in `metricSpec`, so a residual is
expected rather than a bug.

**Version-dependent Bombcell API.** Six spots in this script depend on the installed
Bombcell version and are marked in place: `param.splitGoodAndMua_NonSomatic`,
`param.plotGlobal` / `param.plotDetails`, the `bc.qm.getQualityUnitType` signature, the
RPV classification threshold field (`local_relaxRPV`), the numeric→string label map
(`local_mapNumericBC`), and the per-metric param field names (`local_bcThresholds`).
The header carries an index of all six.

**Runtime:** `forceRerun = false` reloads cached per-recording `bombcell_output/` data
instead of recomputing, so full reruns are fast. The pooled MAT and the per-recording
TSVs must come from **one single full-set run**, so post-promotion labels are consistent
across recordings.

---

## Stage 2 — Unified cluster QC

### `HKS_cluster_QC_unifiedV7.m`

The `_layers` suffix was removed from the filename; older documentation calls this
`HKS_cluster_QC_unifiedV7_layers.m`.

**In:** `kilosort4/` via `loadKSdir`; `cluster_bc_unitType.tsv`; `Rec*ProbeLocation*.mat`
(UniversalProbeFinder — **the `Rec` prefix is required**, a bare `*ProbeLocation*.mat`
matches a deprecated file); `*_stimfile_new.mat`; `kilosort4/syncSY.mat`; NI `.meta`.

**Out:** `<recName>_cluster_QC.mat` — despite the name this is the complete
per-recording preprocessing artifact, not just a mask:

- `sQC.clusterTable` — one row per cluster. IDs and metrics (`clusterID`, `globalID`,
  `nSpikes`, `depth_um`, `probeDepth_um`, `ksGroup`, `bcType`, `wfDuration_us`,
  `cellType`, `firstSpikeLat_ms`, `zetaP_pulse`, `lowbandExcess`, `responseClass`),
  anatomy (`region`, `regionAllen`, `v1Layer`), and per-criterion booleans
  `pass_ksLabel`, `pass_bcLabel`, `pass_notArtefact`, `pass_responsive`.
- `sQC.spikes` — long-format spike list for the whole recording (`spikeTime` on the AP
  clock, `clusterID`), joined to the table on `clusterID`.
- `sQC.stim` — pulse and flicker onsets on the AP clock, so downstream scripts never
  reparse `structEP`.
- `sQC.meta` — every threshold and choice used, including the sync corrections and
  `NS_CUTOFF_US`.

Because the spikes and stimulus live here, downstream scripts load this one file instead
of redoing NI+AP sync, `loadKSdir`, depth extraction and stimfile parsing. That
consolidation is deliberate: it made three historical loading inconsistencies
structurally impossible.

**There is no authoritative include column.** `includeDefault` is a convenience AND of
the pass flags; inclusion is decided at load time by downstream toggles. The operative
analysis gate is **Bombcell-good + pulse-ZETA responsive, KS label off**.

⚠ **`useKS` is off and should stay off.** The Kilosort `group` label proved unreliable
against Bombcell's classification, so `pass_bcLabel` is the authoritative label gate.
`pass_ksLabel` is still stored per cluster and the toggle still exists in the pooler and
both descriptor scripts — it is a documented off switch, not dead code. Turning it on
further restricts the set and will not reproduce the reported numbers. The loader example
at the foot of the QC script now carries these settings.

**Gating decisions worth knowing**

- The **flat pulse-latency artefact screen is an audit column, not a gate**
  (`cfg.INCLUDE_ARTEFACT_GATE = false`). It is rate-confounded — flagged units fire ~3×
  faster, which is an order-statistic effect, not an artefact — it caught ~0 real
  artefacts, and it discarded ~140 good units per recording. A latenZy competitive check
  was built to test this and confirmed it; latenZy has since been removed from
  production QC, its diagnostic job done.
- **Response classification is lowband-only, driven vs suppressed.** Per cluster,
  spontaneous-corrected spikes/cycle averaged over the 4–16 Hz anchor gives
  `lowbandExcess`; sign against floor 0 gives `responseClass = +1` / `−1`. The tuning
  curve spans the full 2–243 Hz; the anchor is only the readout window.
- **Anatomy comes from a `cluster_id` join** on the ProbeLocation file — no depth bands,
  no tip/surface flips. `v1Layer` is parsed from the verbatim Allen name
  (`"…layer 2/3"` → `"L2/3"`) and stored at fine granularity (L1 / L2-3 / L4 / L5 / L6a /
  L6b); grouping is a plot-time decision.
- **`cellType`** is derived at load time from `wfDuration_us` at
  `cfg.NS_CUTOFF_US = 400 µs` (`< 400` narrow-spiking, `≥ 400` broad-spiking), matching
  Bombcell's cortical threshold. Cell-type splits are **V1-only by design** — the cutoff
  is cortex-calibrated.

⚠ **Two things to know before running.** There is no skip-if-exists guard — every
recording in `recSpec` is fully recomputed and overwritten, ZETA included. To re-run one
recording, narrow `recSpec` rather than letting the batch loop. And **ZETA is unseeded**
(250 resamples/cluster), so `zetaP_pulse` shifts slightly between runs and borderline
`pass_responsive` cells can flip, moving the included population by a cell or two. This
was left as-is deliberately; add `rng(7,'twister')` near the top if you want
bit-reproducibility.

**`recSpec` date filtering** is honoured here and in the pooler and the Bombcell batch.
Narrowing `recSpec` re-runs only the named dates; a date with no matching folder warns.

---

## Stage 3 — Stage-1 reduction (pooling)

### `HKS_pool_resultantsV3.m`

The architectural centrepiece. Every following metric is a function of three numbers per
(cluster, frequency), so spikes are only raw material for an `nClust × nFreq` reduction:

```
Zc = mean_k exp(-i·2π·f·t_k)    complex mean resultant   (|Zc| = VS)
Nc = spike count
Tc = summed observation duration (s)
```

plus the matched blank-null triple `Zn / Nn / Tn`. Peak memory stays bounded to one
recording's spikes at a time.

**In:** per-recording `<recName>_cluster_QC.mat`, resolved from `recSpec`.
**Out:** `W:\Kersten\Pooled\resonance_resultants_pool.mat` — `clusterTable` (all kept QC
rows, every column preserved), `freqList` (union grid, **NaN where a recording lacks a
frequency**), `Zc/Nc/Tc`, `Zn/Nn/Tn`, `matMuPeriod`, `vecSpontRate`, `flickerZetaP`
(continuous, not thresholded), `meta`. Also caches a per-recording `_resultants.mat`.

`Zc` is stored as the **mean** resultant, so `Zc.*Nc` recovers the raw complex sum and
both per-cluster averaging and future spike-level re-aggregation stay possible.

**The inclusion gate is baked in here** (`useKS = false`, `useBC = true`,
`useArtefact = false`, `useResponsive = true`). Downstream can only further restrict, so
this is a ceiling, not a floor. `mouseID` and `recID` are retained so animal-weighted
re-pooling stays a downstream choice.

**Analysis windows — the most consequential settings in the pipeline**

- **Real window = the last 3.0 s of the stimulus** (`onset + dblStimDur − dblRealDur`).
  Steady state; the onset transient is excluded.
- **Null window = the 3.0 s immediately before onset.** Durations are matched, which is
  required: unequal real/null durations make the paired `signrank` tests liberal by
  construction, via the count dependence of the PPC estimator.
- **F0 (`matMuPeriod`) uses the same 3.0 s window**, so F0 and |F1| describe the same
  epoch. They sit in adjacent panels of Figure 3, and with adaptation the full-trial mean
  rate exceeds the steady-state tail.
- **Trial selection is symmetric** — a trial enters both sides or neither, so `Nc/Tc` and
  `Nn/Tn` always cover the same trials.
- `Zc` phase is referenced to **stimulus onset**, not window start. A constant offset
  leaves `|Zc|` (and PPC, VS, |F1|) unchanged but rotates `angle(Zc)`, and phase matters
  for any future f/2 cycle-skipping analysis.
- `dblStimDur` is **measured per recording and asserted, never defaulted**, because it
  sets the window *start*: a wrong value slides the analysis window rather than
  mis-scaling a denominator.

**Window headroom.** Flicker trials run ~4.99 s (`meta.dblStimDurPerRec`), separated by
~4 s blanks. Both 3.0 s windows fit inside their epochs with roughly 1 s of slack against
the shortest observed blank gap (4.02 s — measured before M9 was dropped, so **recompute
on the 5-mouse set before relying on it**). Lengthening either window past the shortest
gap trips an assert in the pooler.

⚠ **Two consequences.** `Nc` drops roughly 40% relative to the old full-trial window —
watch guard survival, especially in CA1. And any downstream F0 conversion must divide
`matMuPeriod` by **`meta.dblRealDur` (3.0 s)**, not `meta.dblStimDurPerRec` (4.99 s).
The two thesis metrics scripts are patched for this; anything else reading `matMuPeriod`
has not been checked.

⚠ The pooled-PSTH FFT is **not** reducible to resultants and is not stored. Those
sections would need spikes if ever revived.

---

## Stage 4 — Figures and tables

All read the pooled resultant file (except the rasters). They produce **individual graphs
for hand-stitching** into the composite deck (`composite_figures_updated.pptx`), saved
flat to `W:\Kersten\Figures\Thesis` with the mode and metric in the filename.

Pooling for the descriptive curves is flat cluster mean ± SEM. That SEM is for
**visualisation only** — the inferential unit is the animal-level bootstrap.

### `HKS_thesis_metrics_pooledV3.m` — Figure 3
Three per-area figures: **F0** (spk/s, ungated), **Δ|F1|** (spk/s, null-subtracted,
ungated), **PPC** (Vinck 2010, ungated, spike-guarded). CA1 is kept in all three for a
consistent grid. Output `POOLED_{F0,F1,PPC}_perarea_linkedV1CA1.png`.

**V1 and CA1 share a y-axis within each metric row; dLGN keeps its own.** Autoscaling
every panel independently inflated CA1's noise to fill its frame, which argued against
the claim that CA1 shows an absence. Putting CA1 on the dLGN scale instead flattens V1
into an unreadable line — dLGN is ~5× larger in both PPC and |F1| — so only the two panels
within a few-fold of each other are linked. Every panel keeps its own **tick** labels,
which is what makes the split unambiguous.

**A grey blank-window band sits behind each curve**, in one of three modes, and the
caption must match which:

| mode | metric | what the band is |
|---|---|---|
| `direct` | PPC | the blank-window metric itself, mean ± SEM, in the metric's own units |
| `diff` | Δ\|F1\| | the plotted quantity is already a difference, so there is no second null window. The band is centred on zero and widened by √2 — a **conservative reference scale, not a formal null**, since shared clusters give positive covariance |
| `none` | F0 | already spontaneous-corrected in the reduction, and the pool doesn't carry the uncorrected blank rate |

**Panel furniture.** One y-label per row (`ylabelPanel`), since tick labels already mark
the scale split. Titles are "area (n = X)" only — the metric name lives in the y-label.
The legend sits on `legendPanel` (CA1, the flattest curve) and names **only** the grey
band; the coloured curve is identified by the title. `legendPanel = 0` drops it entirely.

Blank-window PPC is computed in `[1]` rather than in `[N]`, so the band on the figure and
the numbers in the diagnostic table cannot drift apart.

⚠ Where the blank |F1| is undefined it is treated as **zero** rather than masked, so an
|F1| point survives where the matching PPC point would be dropped. Defensible (|F1| is
unguarded by design) but it is an asymmetry between two panels presented as comparable.

**Also carries the blank-null PPC audit** (block `[N]`, output
`POOLED_PPCnull_perarea.png`). PPC is count-debiased, so `E[PPC] = 0` under a uniform
null and the blank curve should sit on zero. It does — which is what licenses plotting raw
PPC rather than ΔPPC everywhere downstream. If that decision is ever revisited, this is
the block to re-run: the `sd ratio` column quantifies the variance cost of subtracting,
and `r(R,N)` explains why the cost is smaller than it looks.

### `HKS_thesis_metrics_supplementV1.m` — Figure 3 supplement
F0 and |F1| in the two non-rate units → `POOLED_{F0,F1}_{percycle,perperiod}_perarea.png`.

**Why the supplement exists.** PPC is dimensionless, so its frequency profile is
unit-invariant — one shape, no freedom. F0 and |F1| are the only rate-carrying metrics,
so the per-second / per-cycle / per-period choice exists only for them, and **|F1|'s
shape is unit-dependent**: dividing spk/s by frequency imposes a monotone 1/f tilt that
turns a band-pass rate curve into an apparent low-pass. dLGN |F1| is band-pass in spk/s
and low-pass in spikes/cycle — the same spikes, opposite qualitative shape. Main figure
carries spk/s, supplement carries per-cycle. Never write "|F1| is low-pass, PPC is
band-pass" without naming the |F1| unit; in spk/s they agree.

Conversions, with `W = meta.dblRealDur = 3.0 s` (**not** the ~4.99 s trial):

| quantity | formula | unit |
|---|---|---|
| F0_period (native) | `matMuPeriod` | spikes / window |
| F0_rate | `matMuPeriod ./ W` | spikes / s |
| F0_cycle | `matMuPeriod ./ (W .* f)` | spikes / cycle |
| F1_rate (native) | `2·Nc·\|Zc\| ./ Tc` | spikes / s |
| F1_cycle | `F1_rate ./ f` | modulated spikes / cycle |
| F1_period | `F1_rate .* W` | modulated spikes / window |

Order is immaterial for the null subtraction: real and null share each factor, so
converting the already-subtracted amplitude equals subtracting the converted ones.

### `HKS_thesis_descriptors_PPCV2.m` — Figure 4
**This is stage 4b.** It does *not* wrap `HKS_following_descriptorsV6`; it supersedes it
(see "Superseded scripts" below). Reduces each group's curve to three descriptors:

- **`f_max`** — peak of the curve.
- **`f_half`** — highest frequency on the descending flank at `halfFrac·peak`
  (amplitude 0.5, or 1/√2 for −3 dB power). Drawn as an L-bracket.
- **`f_cutoff`** — last frequency with significant following: one-sided `signrank`
  (median > 0) per frequency, BH-FDR corrected across the grid, plus a "staying-lost"
  run of 2 non-significant bins.

**Two curves, guarded differently on purpose.** `Mcurve` is **raw PPC** guarded on the
real side only, and drives the plotted curve and the geometric descriptors `f_max` /
`f_half`. `Dprim` is the paired real-minus-blank difference guarded on **both** sides,
and drives only the signrank test behind `f_cutoff`. Blank PPC is statistically zero, so
subtracting it buys no bias correction and only adds spread, worst where the signal is
weakest — but a per-cluster empirical null is still more defensible than testing against
an analytic zero, which is why the paired test keeps it.

Loops `groupMode` over **`area`, `v1layer`, `v1celltype`** (the fourth mode
`v1layertype` is skipped — too thin for animal-level CIs) and plots **PPC only**.
`rng(7)` resets per mode. Output
`following_descriptors_PPC_{area,v1layer,v1celltype}_animal.png`.

**Inference.** `bootUnit = 'animal'`, **n = 5**, `nBoot = 1000`. Resample animals with
replacement, take all of a resampled animal's cells, rebuild groups, re-extract all three
descriptors; between-group contrasts are computed **within each iteration** from the same
resampled set, so they're paired. Cells are nested in probes in animals, so treating cells
as independent inflates effective N (Saravanan, Berman & Sober 2020 is the reference).
No p-values are headlined — the deliverable is the paired CI and whether it excludes zero.
The cell-level `signrank` behind `f_cutoff` is anti-conservative under nesting; say so once.

⚠ `cutoffBootMethod = 'signrank'` (the default, and what the thesis uses) re-runs the
cell-level test on possibly-duplicated resampled cells, giving optimistic variance.
`'curveFloor'` avoids that but is a **different estimator** — don't mix them in one table.

⚠ **`[2b]` guard-survival audit.** `nValid` = paired-valid cells per group × frequency,
`nValidCurv` = real-valid cells behind the plotted curve. NS and BS differ sharply in
rate, so the spike guard can shrink the slower group's high-frequency sample. **An
f_cutoff gap that tracks `nValid` rather than the curve is a selection artefact, not
fidelity.** Check this whenever you report a group difference.

**CA1 is a null result, not a small effect.** Its `f_cutoff` is undefined — not
right-censored. The null gate returns **n.d.** for all three descriptors rather than a
spurious finite value: `max()` of a noise curve always returns *some* frequency, and a
bare 2.0 Hz would read as real low-frequency following. Bootstrap CIs are gated on the
same rule, so a gated group returns no interval and no contrast against it is drawn.

### `HKS_spike_count_figureV2.m` — Figure 2 A–D
Consumes only `Nc`, `Nn` and `Tc`. Four panels: **A** per-cluster spike-count heatmap
(full height, left), **B** yield histogram (carries the figure's single legend),
**C** denominator per frequency on a log-Hz axis, **D** spikes/cycle exposure envelope.

**C is the counting-floor defence** — a guard can only sculpt `f_cutoff` if it thins
cells preferentially near that frequency, so a frequency-flat denominator means the
cutoff reflects genuine fidelity loss rather than attrition. **D's 1/f bound constrains
|F1| but explicitly not PPC**, which pools phase over all cycles and is untroubled by
cells that fire on a minority of them.

**B also shows CA1 yield is not the reason CA1 is null** — CA1 clusters sit in the same
total-spike range as dLGN and V1, so the null cannot be dismissed as insufficient power.

⚠ **`Nn` is the binding side of the paired guard.** Under the 3.0 s null window the
blank-side count drops roughly 40% of bins at high frequency. The paired guard, not the
real-only guard, is what the descriptors apply, so it is what panel C shows.
`showPairedGuard = false` exists only for inspecting the real-trial guard in isolation.

⚠ Two panel-A caveats: it draws **two kinds of missing differently** (NaN `Nc` = frequency
not presented → white; below-guard = grey overlay), and it is on an index axis (`imagesc`
can't do non-uniform spacing) while C and D are on log-Hz — don't compare x-positions
across them.

⚠ This script writes to `W:\Kersten\Figures\` rather than `...\Figures\Thesis\` like the
other stage-4 figures. Change `figDir` if it should sit with them.

### `HKS_unit_rastersV3.m` — Figure 2 E–G
Rasters for hand-picked units, one per region. **The one stage-4 script that does not read
the pooled file** — it loads `<recName>_cluster_QC.mat` per target unit for spikes and
stim times, which guarantees clock consistency with the PPC computation and needs no
NI/CatGT re-sync. Requires `sQC.stim`; it does not reparse `structEP`.

An unresolvable target is **skipped with a warning** rather than aborting, so the
remaining panels still render while a target is being settled.

The deliverable is the **combined composite panel**. Per-unit figures and the short
500 ms window are behind toggles (`cfg.doIndividual`, `cfg.windows`) rather than deleted,
so either comes back without editing `[P1]`. Panel titles are one line — region, mouse,
cluster — since depth and the BC label are QC provenance rather than figure content.

⚠ The frequency axis is capped by `cfg.dblFreqHigh`, currently **40 Hz** — a legibility
choice (above ~50 Hz the bands are too thin to read), not an analysis one. dLGN PPC curves
extend to ~62 Hz, so if a raster and a PPC curve sit side by side that mismatch is visible.

**Exemplars are picked automatically, not by hand.** `HKS_pick_exemplarsV1.m` ranks units
by band-mean ΔPPC within per-area bands (dLGN 10–50 Hz, V1 and CA1 2–30 Hz) under a
stricter-than-group spike guard (`guardEx = 50`, both real and null), gates on band
coverage and breadth, and breaks ties on stop-band contrast. **CA1 is selected as
*typical* rather than top-ranked** — there is no following in CA1 to rank on, so ranking
it like the others would select on the result the panel illustrates. A same-recording
triple finder puts V1 and CA1 in the same session, so the presence-versus-absence contrast
carries no session confound.

⚠ **The picker is not on the Drive** — check the office PC. The three targets are
hardcoded in `cfg.targets`, so the rasters reproduce without it, but the selection rule
does not survive otherwise. Note also that it ran against `HKS_pool_resultantsV2`, so the
exemplars were chosen from the **old analysis window** using ΔPPC and `guardEx = 50`.
Immaterial for illustrative panels, but re-running it against the current pool would not
necessarily return the same three units.

### `HKS_qc_supplement_tableV1.m`, `HKS_sample_tableV1.m`
Descriptive tables for the thesis (unit counts per region × mouse, per V1 layer × mouse,
QC funnel). Counts and proportions only, so no rerun is needed for metric changes — only
if the inclusion gate or cluster set changes. These replaced the earlier
`HKS_results_tablesV1.m`, which is no longer on the drive.

⚠ A responsiveness-proportion table is meaningless on a pool built with
`useResponsive = true` (every unit passes by construction). It needs a pool built with
`useResponsive = false`.

### `HKS_guard_sweep_tablesV1.m` — supplement
Spike-guard sensitivity sweep, `minSpkPerBin = [5, 10, 20, 50]` × three grouping modes,
reference guard **10**. Tables only, no figures. Inclusion, bootstrap unit, estimator and
seed are frozen so only the guard moves; `rng` is re-seeded per (mode, guard) so
differences reflect the guard, not resample noise.

It matches `HKS_thesis_descriptors_PPCV2` — same `Mcurve`/`Dprim` split, same null gate —
which is what makes its numbers comparable to the body. Four outputs:

- `T_desc` — point estimate + CI per (mode, group, guard), with `*_ungated` columns so the
  null gate is inspectable rather than silent
- `T_contrast` — pairwise differences with excludes-zero verdicts
- `T_surv` — paired-valid cell count per frequency
- `T_stab` — spread of the point estimate across the guard set, and deviation from the
  reference guard

**`T_stab.range_over_guards` is what the methods insensitivity claim rests on.** The
console additionally reports whether any pairwise contrast **flips** its excludes-zero
verdict across the sweep — no flips is the strong form of the claim, and is what should be
quoted. A stable range with an unstable verdict would not support the methods line.

**One global guard across all areas is essential**; per-area guards would confound
methodology with biology in the between-area contrasts. The guard is set by hand and the
sweep *licenses* that choice — deriving it from the data would reintroduce the circularity
the sweep exists to remove.

**Result, for the write-up:** dLGN descriptors are invariant across the whole range;
V1 `f_cutoff` holds at 41 Hz [37–52] at every guard; CA1 is degenerate at all guards and
worsens with higher ones — the correct signature of a null.

---

## Metric definitions

All read-outs are functions of `Zc`, `Nc`, `Tc` and the blank-null `Zn`, `Nn`, `Tn`:

- **VS** `= |Zc|`. Upward-biased at low spike counts (`E[VS²] = 1/N` under the null) —
  gate before quoting.
- **PPC** (Vinck 2010) `= (Nc·VS² − 1) / max(Nc−1, 1)`. **Preferred.** Null-centred and
  count-debiased; negative values are noise, **not** anti-locking. The squared scale
  matters when comparing magnitudes: a 5× PPC difference is ~2.4× in resultant length.
- **Preferred phase** `= angle(Zc)`. With `exp(−i·2π·f·t)` this is the phase **lag** of
  the firing peak relative to cycle onset.
- **|F1|** `= 2·Nc·|Zc| / Tc` (spikes/s). Rate-confounded — it has a built-in ~1/f
  envelope, bounded by spikes/cycle = rate/frequency, which makes raw |F1| look low-pass
  even where phase-locking fidelity is band-pass. PPC strips rate out. Null subtraction
  uses the matched triple `2·Nn·|Zn|/Tn`.
- **F0** `matMuPeriod` = spontaneous-corrected spikes per **analysis window** (3.0 s).
- **Spike guard** `minSpkPerBin = 10` on the normalized metrics, applied as a NaN-mask
  before averaging. This is a **variance control, not a bias control** — it is defensible
  because it is applied identically to every group and because the descriptors are
  invariant across the `[5, 10, 20, 50]` sweep.

**PPC gating.** PPC is **ungated (region mask only) plus spike-guarded**. The driven gate
was a VS necessity — VS is upward-biased at low counts. PPC is already count-debiased, so
the gate would only change the estimand, by conditioning inclusion on the low-band *rate*
response. The driven gate stays on VS only.

**Figure 3 PPC is raw (descriptive); Figure 4 uses raw PPC for the curve and ΔPPC for
significance only.** The same metric appears in two roles. Raw PPC drives the plotted
curves and the `f_max` / `f_half` descriptors; ΔPPC drives the `f_cutoff` test. A
null-subtraction audit (block `[N]` of the Figure 3 script) confirmed the blank PPC is
statistically zero, which is what licenses the raw curve.

---

## Shared helpers

Some of these appear byte-identically in several scripts rather than living in a shared
file, so each script runs standalone with no path dependency. **Where that holds, the
invariant is that every copy stays identical — edit all of them or none.** Each such copy
carries a banner line saying so.

### Still identical everywhere

| Helper | Appears in |
|---|---|
| `computeSpontRate` | `cluster_QC_unifiedV7`, `pool_resultantsV3`, `fourierV4`, `fourierV5`, `resfigV8`, `raw_spikes_plotsV7`, `fourier_pooledV1` |
| `meanSemAcrossClusters` | `thesis_descriptors_PPCV2`, `thesis_metrics_pooledV3`, `thesis_metrics_supplementV1`, `fourierV4`, `fourier_pooledV1`, `following_descriptorsV6` |
| `extractCurveDescriptors`, `cutoffCrossing`, `crossingByFloor`, `bhfdr`, `tern` | `thesis_descriptors_PPCV2`, `guard_sweep_tablesV1`, `following_descriptorsV6` |
| `exportFigure` | `thesis_descriptors_PPCV2`, `thesis_metrics_pooledV3` |
| `plotShadedBand` | `thesis_descriptors_PPCV2`, `thesis_metrics_pooledV3`, `thesis_metrics_supplementV1`, and (as the older 4-argument form) `fourierV4`, `fourierV5`, `fourier_pooledV1`, `raw_spikes_plotsV7`, `following_descriptorsV6` |

### ⚠ Forked

**`perAreaFigure` is forked** — nine arguments in `thesis_metrics_supplementV1`, ten in
`thesis_metrics_pooledV3`. The extra `opt` struct carries the null-band overlay and its
three modes, per-axis extent tracking for the linked y-limits, `fs`-driven typography, the
`ylabelPanel` / `legendPanel` furniture, and `exportFigure` in place of `saveas`. The
pooled body is 41 lines to the supplement's 18.

This fork is **intended, not a defect to patch during handover**. Porting it means
bringing the whole `opt` assembly across plus two helpers the supplement doesn't have
(`applyLinkedYLim`, `exportFigure`), and then regenerating and re-checking the supplement
figures. A nine-argument call to the pooled version errors immediately on `opt.fs`, so
there is no silent-failure risk.

One honest cost of leaving it: the main Figure 3 panels have print-sized type and 300 dpi
export and the supplement panels do not, so side by side in the thesis they will not
match. Porting `opt` across is the obvious first improvement if that supplement is ever
regenerated.

Not part of the set, despite appearances:

- **`plotShadedBandLU`** (`spike_count_figureV2`) is a *sibling*, not a copy — a
  median/IQR band is asymmetric and cannot be passed as a single SE vector.
- **`ternary`** (`bombcell_batchV3`) is a separate function from `tern`; they were never
  in step.

Notes on two of them:

- **`computeSpontRate`** — per-cluster rate from the pre-trial blank window
  `[onset − 3.0 s, onset)`, pooled across trials. A sensitivity check at
  W = 0.25/0.5/1/2/3 s showed the median rate stable across a 12× window change
  (3.889 → 3.981 spk/s).
- **`plotShadedBand`** — shaded ±SEM, not `errorbar` (supervisor preference). Drop NaN
  vertices before `fill` (unlike `errorbar`, `fill` can't handle them); patch first, mean
  line on top; set `XScale` **after** plotting. Signature is
  `(x, mu, se, col, alpha, lw, name)`, with the last three optional: a four-argument call
  gives the original `0.2` / `1.3` appearance. An empty `name` keeps the mean line out of
  the legend, which is what stops `legend()` inventing a `data1` placeholder. The
  superseded diagnostics still carry the older four-argument body; that is compatible for
  every call they make, but bring them to this signature if you touch them.

**Figure styling is centralised.** `thesis_descriptors_PPCV2` and
`thesis_metrics_pooledV3` each carry an `fs` struct (axis / label / title / legend sizes
and line widths) plus `figPos` and `exportDPI`. Type is sized for panels that get reduced
into the stitched composite, so it looks oversized on screen on purpose. Export goes
through the shared `exportFigure` helper, which wraps `exportgraphics` and falls back
quietly to `saveas` — `saveas` alone rasterises at screen resolution, which is why labels
vanished once panels were scaled down.

---

## Superseded scripts

### `HKS_following_descriptorsV6.m` — superseded by `HKS_thesis_descriptors_PPCV2.m`

**Do not run it to reproduce reported values.** It is the previous generation of the
descriptor analysis and produces **different numbers from the same pool**:

| | `following_descriptorsV6` | `thesis_descriptors_PPCV2` |
|---|---|---|
| Plotted curve | ΔPPC (null-subtracted) | raw PPC |
| `f_max` / `f_half` from | ΔPPC | raw PPC |
| Significance gate | **absent** | present → n.d. |
| `minSpkPerBin` | 50 | 10 |
| Modes | one per run (4 available) | loops 3 in one pass |

Run V6 on `area` and CA1 returns a finite `f_max` and `f_half` instead of n.d., which
contradicts the reported null. The three changes were each independently reasoned: blank
PPC is statistically zero (so subtracting adds spread), max-of-noise is meaningless (so
gate it), and descriptors are guard-invariant (so the stricter guard costs cells for
nothing).

Its one unique output is the companion **Δ|F1| panel** — the amplitude view alongside
locking fidelity, which is the interpretive payoff of the whole analysis. If that panel is
wanted again, **port it forward rather than running V6**, and keep it **null-subtracted**:
the blank carries real spectral content for |F1| even though blank PPC is zero, so the
argument that moved PPC to a raw curve does not transfer.

### `HKS_pool_resultantsV2.m`
Uses the old analysis window. Running stage 3 with it will not reproduce the reported
numbers. Worth moving to an `outdated/` folder rather than leaving it in a browsable one.

### `HKS_bc_metric_distributions.m`
Everything it plotted is reproduced in BLOCK 5 (panels A and C) of the Bombcell batch, on
the in-memory pooled table. Keep it only for re-inspecting a saved
`bombcell_pooled_QC.mat` without re-running the heavy Bombcell loop.

---

## Diagnostics and side scripts

Off the critical path; none of them feed the main results.

| Script | What it does | Status |
|---|---|---|
| `HKS_raw_spikes_plotsV7` | Per-cluster spike count / peak rate / spikes-per-cycle vs frequency, heatmaps, split-half reproducibility. | Superseded by the pooled path but still useful per recording |
| `HKS_resfigV8` | Per-area FFT resonance spectra + rasters, ProbeLocation-first. | Superseded by the pooled path |
| `HKS_zeta_based_rastersV6` | Flicker ZETA on all clusters, rasters for the top-5 responders, writes `<recName>_zeta_results.mat`. | Working |
| `HKS_fourierV4` / `V5` / `HKS_fourier_pooledV1` / `V2` | The pre-thesis figure path, per recording. | Superseded by the thesis metrics scripts |
| `HKS_pick_exemplarsV1` | Automated exemplar selection for the raster panels (see Figure 2 E–G). Includes a ΔPPC-vs-frequency preview figure and a same-recording triple finder. | **Not on the Drive — recover from the office PC** |
| `HKS_presence_flickerV1` | Recomputes presence restricted to the flicker window, since Bombcell's `presenceRatio` spans the whole recording and penalises units that drift off *after* the flicker block. Flags recoverable units. | Working, not applied |
| `HKS_spatial_footprint_sweepV2` | Template footprint (channel spread, decay) to test whether high-frequency followers are compact (real) or spread (artefact). Frame-independent. | **Complete for M7+M8 — clean** |
| `HKS_post_QCfourier_artefactcheck` | Crosstab of latency suspects × high-frequency following. Runs on the live workspace after `fourierV4` with all mask toggles off. | Done; part of the evidence retiring the latency gate |
| `HKS_rpv_ppc_checkV1` | Tests whether RPV contamination attenuates PPC here. Needs a diagnostic pool with `useKS = useBC = false`, or high-RPV units are already gone. | **Not run** |
| `HKS_ccg_dedup` | Cross-correlogram duplicate / oversplit detector. | **Not run** |

⚠ `HKS_fourierV4` was the origin of the resultant loop that `pool_resultants` copied, but
the V3 window change means the two are **no longer identical** — only the shared helpers
still are.

**Inherited (do not modify):** `RMS1_stimstructure`, `RMS2_resfig`, `RMS3_iCSD`
(current-source-density, never brought into the HKS pipeline — relevant if you pick up
the LFP side), `RMS4_artefacts`, `getPreProSynthesis` (Acquipix), and
`preAP_structEP_patch` (a one-off repair adding missing `ActOnSecs`/`ActOffSecs` fields
to a stim file; **already applied**, hardcoded to one folder, generalise the glob before
any reuse).

⚠ In `RMS2_resfig`, the per-frequency FFT loops must keep running even though their inner
`figure()` calls are commented out. A cleanup that strips the loop on seeing dead
`figure()` calls will silently break `resonance_spectrum` / `ephys_frex`.

---

## Open items

1. **`rpvPromoteCutoff = 0.15`** is hand-picked, partly on yield, and unvalidated. It sets
   the analysed cohort. Sweep it off the saved `bombcell_pooled_QC.mat` if this work goes
   further. See Stage 1.
2. **`useResponsive` at pool time.** `true` bakes pulse-ZETA into the pool, so downstream
   can't recover pulse-unresponsive cells. A phase-locking-only analysis arguably wants
   `false` — a rate-unresponsive cell can still lock — and an honest responsiveness
   denominator needs that pool. Unresolved.
3. **The blank-gap headroom figure (4.02 s)** predates M9's exclusion. Recompute on the
   5-mouse set before quoting it.
4. **`HKS_pick_exemplarsV1.m` is missing from the Drive.** It produced the raster
   exemplars and exists only on the office PC. Same failure mode as the three stage-4
   scripts — worth one sweep of that machine for anything else uncatalogued.
5. **Two table scripts** (`qc_supplement_tableV1`, `sample_tableV1`) have not had the
   documentation cleanup pass.
6. **`perAreaFigure` is forked** between the two Figure 3 scripts (see Shared helpers).
   Recorded as intended. Port `opt` across if the supplement is ever regenerated, so its
   panels match the main figure's typography and export resolution.
7. **`[N.3]` in `thesis_metrics_pooledV3`** is a guarded sanity check that still uses the
   pre-V3 null definition (`onset − spontWindow − stimDur`), i.e. the window that
   overlapped the preceding stimulus. It no longer crashes, but it validates a window the
   pipeline abandoned.
8. **Dead variables**, reported but not removed: `region = 'Cortex'` in the Bombcell batch
   (assigned, never read — the classifier is steered by
   `param.templateDuration_Ctx_threshold`); `metricPrimary` in the descriptor script;
   `recKey` in both metrics scripts. Also `nanmedian` / `nanmin` in the QC script are
   deprecated in current MATLAB, and `meta.bcGoodTypes = {cfg.bcGoodTypes}` double-wraps
   a cell.
9. **Naming residue:** `fakeStimOnTime` / `fakeStarts` in the pooler are named after the
   retired random-ITI null. They now hold the pre-trial blank window. Rename if you ever
   do a rename pass.
10. **Never run to conclusion:** CCG dedup on the kept set; the RPV-tight robustness check
   on the area-PPC gradient; a GOOD-vs-MUA following comparison; white-noise
   reverse-correlation kernels; the f/2 cycle-skipping discriminator (VS at the drive
   frequency cannot distinguish 1:1 locking from 1:2 skipping — the discriminator is VS
   at `f/2`); the field-proxy choice for the single-unit-vs-field question.

**Retired — do not resurrect:** the `catgt_*_copy` folder convention (a temporary measure
to avoid working in a predecessor's folders; no script resolves, prefers or falls back to
one, and any `*_g0_copy*` derived files still on `W:` are orphaned and safe to delete);
depth-band area assignment and tip/surface flips (anatomy comes from the ProbeLocation
join); three-method response classification (`peaksign`, `signrank`); the flat
pulse-latency artefact gate as an inclusion criterion; latenZy in production QC; the
random-ITI fake-onset null (pre-trial-blank only).
