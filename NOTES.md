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

The sensitivity analysis re-fits the de novo mode with `Kmax` at 12/15/20 and
the coefficient-shrinkage prior at three strengths, then matches each solution
one-to-one against the reference by cosine and reports whether the same
COSMIC-matched signatures, activities and coefficients come back. 

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
depart from the model. 


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
