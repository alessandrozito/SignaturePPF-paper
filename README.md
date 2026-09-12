# Poisson process factorization

Code to reproduce figures, simulation results, and real data analysis results from the paper

> Poisson process factorization for modeling mutational processes along cancer
> genomes (2025). Manuscript under review.

## As a first step, install the SignaturePPF package

The package source is included in this repository.

```
install.packages("SignaturePPF_0.1.0.tar.gz", repos = NULL, type = "source")
library(SignaturePPF)
```

It ships with a tutorial covering the model and a worked example on simulated
data:

```
vignette("getting-started", package = "SignaturePPF")
```

Other packages used here: `tidyverse`, `patchwork`, `GenomicRanges`, `rtracklayer`,
`BSgenome.Hsapiens.UCSC.hg19`, `ggalluvial`, `corrplot`, `RcppHungarian`, `RhpcBLASctl`.

Paths live in `config.R`, which every script sources. Set `SIGNATUREPPF_PAPER` if the
repository is not at `~/SignaturePPF-paper`.

The data and the fitted models ship with the repository, so every figure can be
redrawn without refitting. Two exceptions, both over GitHub's file size limit: the
full MCMC chains for the two applications. Figures S5 and S6 need those samplers
rerun. See [data/README.md](data/README.md) and [output/README.md](output/README.md).

Runtimes and the reasoning behind each analysis are in [NOTES.md](NOTES.md).

## Data preprocessing

* Build the ICGC Breast-AdenoCa cohort (2 kb for the applications, 10 kb for the replication)
  - `R/Preprocess_ICGC_BreastAdenoCA.R`
* Build the 80-cancer cohort (Davies et al. 2017) at 10 kb
  - `R/Preprocess_Breast80.R`
  - `R/Split_Breast80_vcfs.R` (only when rebuilding the per-sample VCFs from the release)
* Donor-level clinical annotation for the ICGC cohort
  - `R/Load_PCAWG_clinical.R`

## Code to reproduce the results of the paper

### Figure 1 - Patterns of mutations along the genome

* Figure 1 panels a and b - Aggregate mutations at the Mb scale, and relationship between mutations at the 2kb scale and signal from H3K9me3
  - `R/Figure1_mutation_landscape.R`

### Simulation studies

* Figure S1, Tables S1 and S2 in the Supplementary material - Simulation I, comparison of PPF and baseline NMF
  - `R/Simulation_main.R`
  - `R/Simulation_functions_main.R`
  - `R/Tables_paper.R`

* Figures S2 and S3, Tables S3 and S4 in the Supplementary material - Simulation II, PPF under model misspecification
  - `R/Simulation_misspec.R`
  - `R/Simulation_functions_misspec.R`
  - `R/Reproduce_figures_Simulation_misspec.R`
  - `R/Tables_paper.R`

### Application 1 - de novo signature extraction

* Files to run the model
  - `R/Application_denovo.R`

* Figure 3, and Figures S5 and S12 in the Supplementary material
  - `R/Reproduce_figures_Application_denovo.R`

* Figure S4 in the Supplementary material - signatures under PPF, CompressiveNMF and SignatureAnalyzer
  - `R/FigureS4_signature_comparison.R`

* Figure 2, and Figure S7 in the Supplementary material - model checking and patient-level deviations
  - `R/Figure2_goodness_of_fit.R`

* Figure S8 and Table S5 in the Supplementary material - sensitivity to `K` and to the priors
  - `R/Application_denovo_sensitivity.R` (fits the scenarios)
  - `R/Reproduce_figures_Application_denovo.R` (draws the figure)
  - `R/Tables_paper.R`

### Application 2 - Fixed signatures analysis

* Files to run the model
  - `R/Application_refit.R`

* Figure 4, and Figure S6 in the Supplementary material
  - `R/Reproduce_figures_Application_refit.R`

### Additional analyses in the Supplementary material

* Figure S9 - stability of the covariate effects under forward selection
  - `R/Application_stability_of_covariates.R`

* Figure S10 - replication of the fixed-signature analysis in the Breast80 cohort
  - `R/Application_replicability_80Breast.R`

* Figure S11 - comparison against `TensorSignatures` on chromatin states
  - `./setup_tensorsig_env.sh` (once, to build the Python 3.7 conda environment)
  - `R/Comparison_TensorSignatures.R`
  - `python/run_tensorsignatures.py`, `run_ts_sweep.sh`

* Copy-number composition of the regions flagged by the goodness-of-fit analysis
  - `R/HighCN_arm_composition.R`

## Helper files

Loaded together by `load_functions()`.

| File | What |
|---|---|
| `R/Utils_functions.R` | cohort to model form, mutation assignment, intensity |
| `R/Plot_functions.R` | figure helpers and palettes |
| `R/Preprocess_functions.R` | binning, bin weights, covariate tracks |
| `R/Application_functions.R` | MAP restarts, checkpointed MCMC, summaries |
| `R/GoodnessOfFit_functions.R` | dispersion, time-rescaling tests, residual maps |
| `R/TensorSignatures_functions.R` | chromatin states, tensor export, comparison |
| `R/Simulation_functions.R` | generative core shared by both simulation studies |

## Figure files

Scripts write plots to `figures/`. Several paper figures are assembled from more
than one file:

| Paper figure | Built from |
|---|---|
| Figure 1 | `Figure1_a_total_mutations_Mb.pdf`, `Figure1_b_2kb_mutations_H3K9me3.pdf` |
| Figure 2 | `Figure3_recontructed_Mbscale_join.pdf`, `Figure3_misspecification.pdf` |
| Figure 3 | `Figure4_a_b_c_denovoPars.pdf`, `Figure4_d_denovoBaselines.pdf` |
| Figure 4 | `Figure5_a_b.pdf`, `Figure5_a_mu.pdf`, `Figure5_c_.png` |
| Figure S1 | `Simualations_results_Supplement2.pdf` |
| Figure S2 | `Simulation_misspec_recovery.pdf` |
| Figure S3 | `Simulation_misspec_calibration_curves_ECE.pdf` |
| Figure S4 | `Breast_suppl_Signatures_comparison.pdf` |
| Figure S5 | `FigureS4_ess_denovo.pdf` |
| Figure S6 | `FigureS5_ess_refit.pdf` |
| Figure S7 | `FigureS_residuals_along_genome.pdf`, `FigureS_residuals_by_mutation_class.pdf`, `GoF_count_distribution_2kb.pdf` |
| Figure S8 | `Sensitivity_mu_all_scenarios.pdf`, `Sensitivity_signatures_comparison.pdf` |
| Figure S9 | `Stability_betas_sequence.pdf`, `Stability_riverplots.pdf`, `Stability_megabase_residuals.pdf` |
| Figure S10 | `Breast80_BreastICGC_descripion.pdf`, `Replication_mu_both_cohorts.pdf`, `Replication_betas_scatter.pdf` |
| Figure S11 | `TensorSignatures_top_panel.pdf`, `TensorSignatures_top_panel_mu.pdf`, `TensorSignatures_chromatin_effects.pdf` |
| Figure S12 | `FigureS2_covariate_correlation.pdf` (the two pages, side by side) |

Tables are written ready to `\input`:

| Paper table | File |
|---|---|
| Table S1 | `output/Simulation_main/table_MAP_timing.tex` |
| Table S2 | `output/Simulation_main/table_MCMC_ess.tex` |
| Table S3 | `output/Simulation_misspec/table_MAP_timing.tex` |
| Table S4 | `output/Simulation_misspec/table_MCMC_ess.tex` |
| Table S5 | `output/Application_denovo_sensitivity/table_sensitivity.tex` |
