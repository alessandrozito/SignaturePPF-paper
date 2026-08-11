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

## Layout

```
config.R                              paths and shared settings, sourced by every script

R/Application_replicability_80Breast.R    the two-cohort replication analysis
R/Application_stability_of_covariates.R   the nested covariate-set analysis
R/Comparison_TensorSignatures.R       the TensorSignatures comparison, end to end

R/Utils_functions.R                   cohort -> model form, mutation assignment, intensity
R/Plot_functions.R                    figure helpers
R/TensorSignatures_functions.R        chromatin states, tensor export, comparison

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
