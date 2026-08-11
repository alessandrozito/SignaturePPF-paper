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

All paths derive from two environment variables, both with sensible defaults, so
the pipeline moves to a cluster without editing any script:

| Variable | Default | What |
|---|---|---|
| `SIGNATUREPPF_PAPER` | `~/SignaturePPF-paper` | this repository |
| `SIGNATUREPPF_DATA` | `~/SigPoisProcess/data` | preprocessed cohort data |
| `SIGNATUREPPF_CHROMHMM` | `~/E028_15_coreMarks_dense.bed` | Roadmap ChromHMM segmentation |
| `TENSORSIG_PYTHON` | `~/miniconda3/envs/tensorsig/bin/python` | TensorSignatures interpreter |

`SIGNATUREPPF_DATA` still points into the old project because the preprocessed
`data` objects were built there. It changes once `SignaturePPF_preprocess()`
exists and those objects are rebuilt into `data/`.

## Analyses

### 1. Replication across two breast cohorts

```
Rscript R/01_replication_Breast80_vs_ICGC.R
```

Refits the same fixed COSMIC signature set to two independent breast cohorts —
80 breast cancers (Davies et al. 2017) and ICGC Breast-AdenoCa — both binned at
10 kb with the same 11 covariates, and asks whether the estimated genomic
covariate effects agree. The signatures are held fixed, so the only thing being
compared is β, plus which signatures the compressive prior keeps.

MAP only. Runtime a few minutes per cohort.

### 2. Comparison against TensorSignatures

Three steps, because the middle one runs under a different Python.

```
Rscript R/02_tensorsignatures_prepare.R    # build data, fit PPF, export tensor
bash/setup_tensorsig_env.sh                # once: build the conda environment
bash/run_ts_sweep.sh                       # fit TensorSignatures over ranks
Rscript R/03_tensorsignatures_compare.R    # import and compare
```

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
interpreter with anything modern. `bash/setup_tensorsig_env.sh` builds it under
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
config.R                    paths and shared settings, sourced by every script
R/functions/                helpers, loaded by load_functions()
  data_adapters.R           cohort -> model form; intensity reconstruction
  assignment.R              per-mutation signature attribution
  plot_helpers.R            figure helpers
  tensorsignatures.R        chromatin states, tensor export, comparison
R/0*.R                      analysis scripts, numbered in run order
python/                     TensorSignatures runner (Python 3.7 environment)
bash/                       environment setup and the rank sweep
data/                       inputs (not tracked)
output/                     results (not tracked)
figures/                    figures (not tracked)
```

`data/`, `output/` and `figures/` are gitignored: the cohort data is
access-controlled and cannot be redistributed, and everything else is
reproducible from the scripts.

## Figure map

| Figure | Produced by |
|---|---|
| `01_burden_along_genome.pdf` | `R/01_replication_Breast80_vs_ICGC.R` |
| `01_betas_Breast80.pdf`, `01_betas_BreastICGC.pdf` | `R/01_replication_Breast80_vs_ICGC.R` |
| `01_beta_replication.pdf`, `01_beta_difference.pdf` | `R/01_replication_Breast80_vs_ICGC.R` |
| `02_ppf_chromatin_betas.pdf` | `R/02_tensorsignatures_prepare.R` |
| `03_ts_rank_sweep.pdf` | `R/03_tensorsignatures_compare.R` |
| `03_cosine_to_cosmic.pdf` | `R/03_tensorsignatures_compare.R` |
| `03_chromatin_effects_*.pdf` | `R/03_tensorsignatures_compare.R` |
| `03_mutation_rate_along_genome.pdf` | `R/03_tensorsignatures_compare.R` |

## Note on the predecessor package

These analyses were originally written against `SigPoisProcess`, which offered
two parameterisations: the original prior and the activity prior. SignaturePPF
implements **only** the activity prior, so results here are not numerically
comparable to fits made with `SigPoisProcess()` and the old
original-vs-activity comparison cannot be rerun from this repository.
