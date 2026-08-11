# SignaturePPF-paper

Reproducibility materials for the SignaturePPF paper: Poisson process
factorization for mutational signature analysis with genomic covariates.

The method itself lives in a separate repository,
[SignaturePPF](https://github.com/alessandrozito/SignaturePPF). This one holds
only the analyses, and pins the exact package version used for the published
results.

## Setup

```r
# install.packages("remotes")
remotes::install_github("alessandrozito/SignaturePPF")
```

Other R packages used here: tidyverse, GenomicRanges, rtracklayer, BSgenome and
BSgenome.Hsapiens.UCSC.hg19, patchwork, ggalluvial, RcppHungarian, RhpcBLASctl.

`config.R` caps BLAS to a single thread. That is not a throttle: the linear
algebra here is tall-and-skinny, so one thread is within ~15% of the best
setting on wall time while using a tenth of the CPU, and the unlimited default
is actually *slower* than one thread on a 24-core machine. Override with
`SIGNATUREPPF_BLAS_THREADS=8` if you are running a single fit and want the last
15%.

All inputs live in `data/` inside this repository, so the analyses depend on
nothing outside it. Two environment variables override the defaults if needed:

| Variable | Default | What |
|---|---|---|
| `SIGNATUREPPF_PAPER` | `~/SignaturePPF-paper` | this repository |
| `SIGNATUREPPF_DATA` | `<repo>/data` | inputs, if they must live on another volume |
| `TENSORSIG_PYTHON` | `~/miniconda3/envs/tensorsig/bin/python` | TensorSignatures interpreter |

The data files are not tracked — the ICGC cohort is access-controlled — but
`data/MANIFEST.tsv` is, and records the size and MD5 of every file the published
results were computed from. See [data/README.md](data/README.md) for what each
one is and where the public ones come from.

## Analyses

### 0. Build the cohort object

```
Rscript R/Preprocess_ICGC_BreastAdenoCA.R          # 2 kb, what the applications use
Rscript R/Preprocess_ICGC_BreastAdenoCA.R 10000    # the coarser grid
```

Bins the genome, computes usable sequence per bin (assembly gaps and the ENCODE
blacklist removed), averages the eleven covariate tracks onto those bins,
winsorises and standardises them, attaches each mutation's covariate values, and
multiplies copy number by usable sequence to give the exposure the Poisson
process integrates over. About ten minutes and 8 GB at 2 kb. Skipped if the
output already exists.

**One fix relative to the predecessor's loader.** It pre-filled the mutation
covariate matrix with zeros and wrote only the rows that matched a retained bin.
The covariates are standardised, so a zero row is not "missing" — it reads as a
perfectly average bin, and mutations in assembly gaps were silently fitted as if
they sat in one. They are now dropped, which is the only consistent choice: the
model integrates its intensity over the retained bins, so a mutation outside
them has no exposure behind it. The consequence is that
`ICGC_BreastAdenoCA_avg10kb_*.rds.gzip` as shipped was built by the old code and
rebuilding it here will not reproduce it byte for byte.

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

The environment is a separate command on purpose: `setup_tensorsig_env.sh`
downloads and installs miniconda under `$HOME`, which an analysis script should
not do behind your back.

Both methods are fitted to the same mutations on the same genomic partition, the
ChromHMM 15-state annotation of breast epithelium (Roadmap E028). The bins *are*
the ChromHMM segments, so the state assignment is exact for both, and PPF is
given the states as one-hot covariates with `Quies` as reference — which makes
each β a log enrichment relative to `Quies`, the same quantity TensorSignatures
reports as a state amplitude.

The comparison then runs at three levels: signature spectra (matched one-to-one
by cosine similarity), chromatin-state effects, and predicted regional mutation
rate in 1 Mb windows. The third is where the models genuinely differ —
TensorSignatures predicts a total per (state, sample) and has no notion of
position within a state.

The chromatin-state effects are shown cut two ways, side by side: one panel per
matched signature, and one panel per chromatin state. PPF signatures with
`mu < 0.05` are dropped *before* the matching — the compressive prior has
switched them off, so their coefficients are prior draws, and because the
matching is one-to-one a dead signature left in the pool can claim a
TensorSignatures signature and displace a live one.

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
per-bin residual. It answers three questions the reviewers raised: do the
coefficients stay put as further covariates enter, does mutation attribution
move between signatures, and does the extra structure predict on held-out
genome or only fit the training bins.

MAP only, and inherently serial — the covariate chosen at step *m*+1 is a
function of the fit at step *m*, so there is nothing to parallelise. Each of the
twelve fits is cached to its own file, so an interrupted run resumes.

The null model is the same PPF with β held at exactly zero, not a separate NMF:
same likelihood, same compressive prior, same optimiser, so the first point of
every curve is the nested null of the ones after it rather than another method's
answer.

Runs on the full cohort — no samples are excluded.

**Do not report TensorSignatures strand-asymmetry results from this pipeline.**
The mutation channels are already pyrimidine-normalised by the preprocessing, so
the strand orientation each mutation had is lost. Both strand axes are filled
with the "unknown" state, which reduces TensorSignatures to its genomic-state
part — exactly the part PPF can also express, so the comparison stays
apples-to-apples, but the strand results are not there to be read. Recovering
them means redoing the annotation from the raw calls.

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

Both applications find a mode and then sample from it: 10000 sweeps, 5000
burn-in, checkpointed every 100 so an interrupted chain resumes bit-exactly.
The de novo one searches from three random starts and keeps the highest log
posterior, since that posterior is multimodal; the refit needs only one, the
signatures being fixed. Neither prunes the signatures the compressive prior
parks near `epsilon` — which of them survive is a result, not a setting.

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

## Layout

```
config.R                              paths and shared settings, sourced by every script

R/Preprocess_ICGC_BreastAdenoCA.R         build the binned cohort from raw tracks
R/Application_denovo.R                    K = 12 estimated de novo, MAP then MCMC
R/Application_refit.R                     15 COSMIC signatures held fixed
R/Application_denovo_sensitivity.R        stability to Kmax and the priors, MAP only
R/Application_replicability_80Breast.R    the two-cohort replication analysis
R/Application_stability_of_covariates.R   the nested covariate-set analysis
R/Comparison_TensorSignatures.R       the TensorSignatures comparison, end to end
R/Simulation_main.R                   the main simulation study (Section 4)
R/Simulation_misspec.R                the misspecification study, three stages

R/Utils_functions.R                   cohort -> model form, mutation assignment, intensity
R/Plot_functions.R                    figure helpers
R/Preprocess_functions.R              binning, bin weights, covariate tracks
R/Application_functions.R             MAP restarts, checkpointed MCMC, summaries
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
other three scripts are the analyses.

`data/`, `output/` and `figures/` are gitignored: the cohort data is
access-controlled and cannot be redistributed, and everything else is
reproducible from the scripts.

## Figure map

| Figure | Produced by |
|---|---|
| `Replication_burden_along_genome.pdf` | `R/Application_replicability_80Breast.R` |
| `Replication_betas_Breast80.pdf`, `Replication_betas_BreastICGC.pdf` | `R/Application_replicability_80Breast.R` |
| `Replication_betas_scatter.pdf`, `Replication_betas_difference.pdf` | `R/Application_replicability_80Breast.R` |
| `Stability_betas_sequence.pdf` | `R/Application_stability_of_covariates.R` |
| `Stability_riverplots.pdf` | `R/Application_stability_of_covariates.R` |
| `Stability_rmse.pdf` | `R/Application_stability_of_covariates.R` |
| `TensorSignatures_PPF_chromatin_betas.pdf` | `R/Comparison_TensorSignatures.R` |
| `TensorSignatures_chromatin_effects.pdf` (the two cuts side by side), `..._by_signature.pdf`, `..._by_state.pdf` | `R/Comparison_TensorSignatures.R` |
| `TensorSignatures_mutation_rate_along_genome.pdf` | `R/Comparison_TensorSignatures.R` |

## Note on the predecessor package

These analyses were originally written against `SigPoisProcess`, which offered
two parameterisations: the original prior and the activity prior. SignaturePPF
implements **only** the activity prior, so results here are not numerically
comparable to fits made with `SigPoisProcess()` and the old
original-vs-activity comparison cannot be rerun from this repository.
