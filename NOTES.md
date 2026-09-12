# Poisson process factorization for modeling mutational processes along the genome

This folder reproduces the paper's results. 

Methodological detail, runtime figures and design decisions behind each script.
The figure-by-figure index lives in [README.md](README.md).

The method itself is a separate R package, whose source is included here as
`SignaturePPF_0.1.0.tar.gz`. This repository holds only the analyses, and pins
the exact package version used for the published results.

## Setup

```r
install.packages("SignaturePPF_0.1.0.tar.gz", repos = NULL, type = "source")
```

Other R packages used here: tidyverse, GenomicRanges, rtracklayer, BSgenome and
BSgenome.Hsapiens.UCSC.hg19, patchwork, ggalluvial, RcppHungarian, RhpcBLASctl.

All inputs live in `data/` inside this repository. Two environment variables override the defaults if needed:

| Variable | Default | What |
|---|---|---|
| `SIGNATUREPPF_PAPER` | `~/SignaturePPF-paper` | this repository |
| `SIGNATUREPPF_DATA` | `<repo>/data` | inputs |
| `TENSORSIG_PYTHON` | `~/miniconda3/envs/tensorsig/bin/python` | TensorSignatures interpreter |

The data files are not tracked because they are large, but all of them are public
downloads; `data/MANIFEST.tsv` is tracked, and records the size and MD5 of every
file the published results were computed from. See [data/README.md](data/README.md) for what each
one is and where the public ones come from.

## Analyses

### 0. Build the cohort object

```
Rscript R/Preprocess_ICGC_BreastAdenoCA.R          # 2 kb, what the applications use
Rscript R/Preprocess_ICGC_BreastAdenoCA.R 10000    # the coarser grid
Rscript R/Preprocess_Breast80.R                    # the 80-cancer cohort, 10 kb
```

Bins the genome, computes usable sequence per bin (assembly gaps and the ENCODE
blacklist removed), averages the eleven covariate tracks onto those bins,
winsorises and standardises them, attaches each mutation's covariate values, and
multiplies copy number by usable sequence to give the exposure the Poisson
process integrates over. Skipped if the output already exists. About ten minutes
and 8 GB for ICGC at 2 kb; two minutes for the 80-cancer cohort.

The 80-cancer cohort reads per-sample VCFs, but the public release ships its
substitutions as one combined table. `Rscript R/Split_Breast80_vcfs.R` turns that
download into `data/SNP80Breast/`; only needed if rebuilding it from source.

The two cohorts differ only at the front end. ICGC arrives as an assembled
`GRanges` plus one consensus copy-number table; the 80-cancer cohort arrives as
80 CaVEMan VCFs and 80 ASCAT segment tables, so its SNVs are read, filtered to
clean single-base substitutions, and assigned a trinucleotide channel from hg19
first. The remaining steps are shared between the two cohorts.

Mutations that fall outside the retained bins are dropped rather than given a
zero covariate row. The covariates are standardised, so a zero row would be read
as an average bin rather than as missing, and the model integrates its intensity
only over the retained bins, so a mutation outside them carries no exposure.

### 1. Replication across two breast cohorts

```
Rscript R/Application_replicability_80Breast.R
```

Refits the same fixed COSMIC signature set to two independent breast cohorts —
80 breast cancers (Davies et al. 2017) and ICGC Breast-AdenoCa — both binned at
10 kb with the same 11 covariates, and asks whether the estimated genomic
covariate effects agree. The signatures are held fixed, so the only thing being
compared is β, plus which signatures the compressive prior keeps.

MAP only. Runtime a few minutes per cohort.

### 2. Comparison against TensorSignatures

```
./setup_tensorsig_env.sh                  # once: build the conda environment
Rscript R/Comparison_TensorSignatures.R   # everything else, end to end
```

The second command builds the shared dataset, fits SignaturePPF, exports the
tensor, drives the TensorSignatures rank sweep (shelling out to the Python 3.7
environment) and runs the comparison. Every expensive step is cached, so a
rerun only redoes what is missing and an interrupted sweep resumes.

The environment is built by a separate command because `setup_tensorsig_env.sh`
downloads and installs miniconda under `$HOME`.

Both methods are fitted to the same mutations on the same genomic partition, the
ChromHMM 15-state annotation of breast epithelium (Roadmap E028). The bins *are*
the ChromHMM segments, so the state assignment is exact for both, and PPF is
given the states as one-hot covariates with `Quies` as reference — which makes
each β a log enrichment relative to `Quies`, the same quantity TensorSignatures
reports as a state amplitude.

The comparison then runs at three levels: signature spectra (matched one-to-one
by cosine similarity), chromatin-state effects, and predicted regional mutation
rate in 1 Mb windows. The models differ most on the third: TensorSignatures
predicts a total per (state, sample) and has no notion of position within a
state.

The chromatin-state effects are shown cut two ways, side by side: one panel per
matched signature, and one panel per chromatin state. PPF signatures with
`mu < 0.05` are dropped *before* the matching — the compressive prior has
switched them off, so their coefficients are prior draws, and because the
matching is one-to-one, a switched-off signature left in the pool can be matched
to a TensorSignatures signature and displace an active one.

**TensorSignatures needs its own environment.** Version 0.5.0 pins
`tensorflow <= 1.15`, whose wheels stop at Python 3.7, so it cannot share an
interpreter with anything modern. `setup_tensorsig_env.sh` builds it under
`$HOME` with no root access; `rm -rf ~/miniconda3` undoes it.

### 3. Stability of the covariate effects

```
Rscript R/Application_stability_of_covariates.R
```

Holds out a fifth of the genome — whole megabases, stratified by chromosome —
and fits a sequence of nested models on the rest: first with no covariate
effect at all, then adding one covariate at a time by forward selection on the
per-bin residual. Three quantities are tracked: whether the coefficients remain
stable as further covariates enter, whether mutation attribution moves between
signatures, and whether the additional structure predicts on held-out genome or
only improves fit on the training bins.

MAP only, and inherently serial — the covariate chosen at step *m*+1 is a
function of the fit at step *m*, so there is nothing to parallelise. Each of the
twelve fits is cached to its own file, so an interrupted run resumes.

The null model is the same PPF with β held at exactly zero, not a separate NMF:
same likelihood, same compressive prior, same optimizer, so the first point of
every curve is the nested null of those after it.

Runs on the full cohort — no samples are excluded.

**Do not report TensorSignatures strand-asymmetry results from this pipeline.**
The mutation channels are already pyrimidine-normalised by the preprocessing, so
the strand orientation each mutation had is lost. Both strand axes are filled
with the "unknown" state, which reduces TensorSignatures to its genomic-state
part, which is the part PPF also expresses. The comparison remains valid, but
the strand results are not present. Recovering them requires redoing the
annotation from the raw calls.

### 4. Simulation study under misspecification

```
Rscript R/Simulation_misspec.R          # generate, fit, score
Rscript R/Simulation_misspec.R fit      # or one stage at a time
```

Seven scenarios, each adding one violation of the model's assumptions on top of
the previous one — noisy patient-specific epigenome, hypermutation hotspots,
noisy copy number, channel-specific opportunity — with 20 replicates each. The
last 4000 of 20000 tiles are held out of every fit, so every metric is reported
in and out of sample. SignaturePPF (MAP and MCMC) is scored against
CompressiveNMF and SignatureAnalyzer on reconstruction, and against
CompressiveNMF on attribution and calibration with the signatures held fixed.

140 jobs run 20 at a time, one core each and nothing nested. Run it pinned:

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  taskset -c 0-19 Rscript R/Simulation_misspec.R
```

Every dataset and every fit is skipped if already on disk, so an interrupted run
resumes. Budget most of a day with the MCMC on; `RUN_MCMC <- FALSE` cuts it to
about an hour and keeps every metric except the credible-interval ones.

### 5. The two main applications, at 2 kb

```
Rscript R/Application_denovo.R              # K = 12 estimated from the data
Rscript R/Application_refit.R               # 15 COSMIC signatures held fixed
Rscript R/Application_denovo_sensitivity.R  # is the de novo solution stable?
```

Both applications find a mode and then sample from it: 10000 iterations, 5000
burn-in, checkpointed every 100 so an interrupted chain resumes bit-exactly.
The de novo one searches from three random starts and keeps the highest log
posterior, since that posterior is multimodal; the refit needs only one, the
signatures being fixed. Neither prunes the signatures the compressive prior
parks near `epsilon`; which of them survive is an output of the fit.

The sensitivity analysis re-fits the de novo mode with `Kmax` at 12/15/20 and
the coefficient-shrinkage prior at three strengths, then matches each solution
one-to-one against the reference by cosine and reports whether the same
COSMIC-matched signatures, activities and coefficients come back. **MAP only** —
the question is whether the mode moves, and answering it does not need a
posterior at a day per scenario.

At 2 kb this is 1.39 M bins and 707 k mutations, so the chains are runs of many
hours. Launch them detached and pinned:

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nohup setsid \
  taskset -c 0 Rscript R/Application_denovo.R \
  > output/Application_denovo/run.log 2>&1 < /dev/null &
```

### 6. Goodness of fit and patient-level deviations

```
Rscript R/Figure2_goodness_of_fit.R
```

Model checking on the 2 kb de novo fit, addressing two questions: whether the
fitted intensity reproduces the observed counts and whether the residuals are
systematic by region, chromosome or mutation type; and where individual patients
depart from the model. **No patient is excluded and nothing is refitted** — the
analysis locates the departures rather than removing the data producing them.

No posterior predictive check. Under the model the counts in disjoint bins are
independent Poisson with a known mean, so everything a predictive simulation
would estimate is available in closed form, including the sampling variability
the check has to be judged against. The distribution of counts at 2 kb, for
instance, is the model-implied marginal
`N_m = sum_bj Pois(m; lambda_bj)`, whose exact variance is
`sum_bj p_m (1 - p_m)` because the cells are independent Bernoulli indicators.

**The diagnostics report effect sizes, not p-values.** With 707k mutations on
1.39M bins every consistent test rejects for deviations far below the size that
would change a conclusion, so a p-value here measures the sample size and not
the model. Three substitutions are made:

- The **dispersion ratio** `D = X²/R` estimates a fixed quantity and does not
  grow with `n`. If the regional rate carries a multiplicative error of
  coefficient of variation `cv`, then
  `D = 1 + cv²·mean(e)`, so `cv = sqrt((D-1)/mean(e))` is the fractional error
  in the predicted rate beyond Poisson noise. `D` is judged against its **exact**
  Poisson sampling band (`Var((X-e)²/e) = 2 + 1/e`) rather than chi-square
  asymptotics, so it stays valid at 2 kb where the expected counts are far below
  1, and no minimum expected count is imposed anywhere.
- Where a test is still wanted, **power is equalised by Poisson thinning**. Each
  patient is binomially thinned to the same expected event count before the
  time-rescaling test; thinning a Poisson process leaves a Poisson process, so
  the test stays exact while a 64k-mutation hypermutator and a 1.2k-mutation
  patient become comparable. Without it the per-patient p-values simply order
  the patients by burden, and the apparent association between hypermutation and
  misfit is an artefact of power alone.
- Rejection is replaced by **localisation**: `X²` is decomposed across regions,
  and what is reported is how much of it sits in the worst 0.1% of them and what
  `D` falls to once they are removed. This separates excess dispersion carried by
  a small number of hotspots from misfit spread across the genome.

The time-rescaling test uses the **conditional** (uniform) form rather than the
exponential-gap form. The total is a fitted quantity the model matches by
construction, so a test that also charges for it spends power on something
already known to agree. What remains is the regional shape, which is the
quantity the covariates model. Within-bin positions are randomised (a randomised
PIT), because the fitted intensity is piecewise constant on a bin and the model
never claimed within-bin structure.

**Two per-patient tests are reported, and they answer different questions.**
`p_full` runs the time-rescaling test on every mutation a patient has, and is
what says whether that patient is misspecified — a patient with 64,000 mutations
genuinely carries more evidence than one with 1,200, and discarding it to make
them comparable discards real power. `p_thin` runs it after thinning to `n0`, and
is *only* for asking whether misfit is associated with burden, where that power
difference is the confounder.

Thinning does not shrink a deviation, it raises the bar: a spike of *k* events
moves the empirical CDF by about `k/n`, which thinning leaves alone, while the
critical distance grows like `1/sqrt(n)`. A focal excess that is decisive at full
`n` can therefore sit under the threshold after thinning — correct for a
power-matched comparison, wrong for detection. On this cohort the two arms give
86 and 36 rejections of 113.

When the thinned arm is used, `n0` must be at or below the smallest burden in the
cohort: a patient with fewer events cannot be thinned up to it, keeps all of its
events, and is then tested at *higher* resolution than the rest while appearing
in the same figure. `n0 = NULL` (the default) sets it to the cohort minimum; an
explicit larger value warns and names how many patients it affects.

Every step caches to `output/GoodnessOfFit/`. About fifteen minutes and 12 GB
from cold; MAP only, so it needs no chain.

`R/Figure2_goodness_of_fit.R` builds the paper figure from the same machinery:
the megabase reconstruction track (with the covariate-free prediction overlaid),
the two reconstruction scatters, the 1 Mb Pearson residuals with each chromosome's
rejection rate written along the top, and a
patient-level volcano (KS distance against significance, on all of each
patient's mutations). The equal-power QQ and the (patient, megabase) cell
volcano appear in the supplement. Each flagged region is also decomposed over
patients, and a region whose excess comes from a single patient is drawn open
rather than filled.

## Layout

```
config.R                              paths and shared settings, sourced by every script

R/Preprocess_ICGC_BreastAdenoCA.R         build the ICGC cohort from raw tracks
R/Preprocess_Breast80.R                   build the 80-cancer cohort from VCFs + ASCAT
R/Split_Breast80_vcfs.R                   split the released CaVEMan table into those VCFs
R/Application_denovo.R                    K = 12 estimated de novo, MAP then MCMC
R/Application_refit.R                     15 COSMIC signatures held fixed
R/Application_denovo_sensitivity.R        stability to Kmax and the priors, MAP only
R/Application_replicability_80Breast.R    the two-cohort replication analysis
R/Application_stability_of_covariates.R   the nested covariate-set analysis
R/Figure2_goodness_of_fit.R            model checking, Figure 2 and its supplements
R/Comparison_TensorSignatures.R       the TensorSignatures comparison, end to end
R/Simulation_main.R                   the main simulation study (Section 4)
R/Simulation_misspec.R                the misspecification study, three stages

R/Utils_functions.R                   cohort -> model form, mutation assignment, intensity
R/Plot_functions.R                    figure helpers
R/Preprocess_functions.R              binning, bin weights, covariate tracks
R/Application_functions.R             MAP restarts, checkpointed MCMC, summaries
R/GoodnessOfFit_functions.R           dispersion, rescaling tests, residual maps
R/TensorSignatures_functions.R        chromatin states, tensor export, comparison
R/Simulation_functions.R              the generative core of the simulation study
R/Simulation_functions_misspec.R      misspecification scenarios, fitting, scoring
R/Simulation_functions_main.R         main-study generation, coarsening, scoring

setup_tensorsig_env.sh                build the TensorSignatures conda environment
run_ts_sweep.sh                       the TensorSignatures rank sweep, driven by the
                                      script above but runnable on its own
python/run_tensorsignatures.py        run under that environment, not the system Python

data/                                 inputs (not tracked; see data/MANIFEST.tsv)
output/                               results (not tracked)
figures/                              figures (not tracked)
```

`R/*_functions.R` are the helpers, loaded together by `load_functions()`; the
remaining scripts are the analyses.

The input data and the fitted models are tracked, so the figures can be redrawn
without refitting. The full MCMC chains, the per-replicate simulation output and
`figures/` are not; see `output/README.md`.
