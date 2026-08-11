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
BSgenome.Hsapiens.UCSC.hg19, patchwork, RcppHungarian.

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

**TensorSignatures needs its own environment.** Version 0.5.0 pins
`tensorflow <= 1.15`, whose wheels stop at Python 3.7, so it cannot share an
interpreter with anything modern. `setup_tensorsig_env.sh` builds it under
`$HOME` with no root access; `rm -rf ~/miniconda3` undoes it.

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

R/Application_replicability_80Breast.R  the two-cohort replication analysis
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
| `TensorSignatures_PPF_chromatin_betas.pdf` | `R/Comparison_TensorSignatures.R` |
| `TensorSignatures_rank_sweep.pdf` | `R/Comparison_TensorSignatures.R` |
| `TensorSignatures_cosine_to_cosmic.pdf` | `R/Comparison_TensorSignatures.R` |
| `TensorSignatures_chromatin_effects_pooled.pdf`, `..._by_signature.pdf` | `R/Comparison_TensorSignatures.R` |
| `TensorSignatures_mutation_rate_along_genome.pdf` | `R/Comparison_TensorSignatures.R` |

## Note on the predecessor package

These analyses were originally written against `SigPoisProcess`, which offered
two parameterisations: the original prior and the activity prior. SignaturePPF
implements **only** the activity prior, so results here are not numerically
comparable to fits made with `SigPoisProcess()` and the old
original-vs-activity comparison cannot be rerun from this repository.
