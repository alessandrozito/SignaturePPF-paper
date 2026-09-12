################################################################################
# Produces: everything. Runs the pipeline end to end, in dependency order.
#
# Usage:
#   Rscript run_all.R                 # every stage
#   Rscript run_all.R preprocess      # one stage, or several
#   Rscript run_all.R figures tables
#   Rscript run_all.R --list          # show the stages and exit
#
# Stages are skipped when their outputs already exist, so a rerun only redoes
# what is missing. The two applications and the two simulation studies run for
# many hours; see NOTES.md for runtimes and for how to launch them detached.
################################################################################

source(if (file.exists("config.R")) "config.R" else "../config.R")

run <- function(script) {
  path <- file.path(R_DIR, script)
  if (!file.exists(path)) stop("no such script: ", path, call. = FALSE)
  message("\n", strrep("=", 78), "\n== ", script, "\n", strrep("=", 78))
  t0 <- Sys.time()
  system2("Rscript", path, stdout = "", stderr = "")
  message(sprintf("-- %s done in %.1f min", script,
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

STAGES <- list(

  ## R/Split_Breast80_vcfs.R is not run here: it rebuilds data/data_for_breast80/
  ## SNP80Breast/ from the released CaVEMan table, and those VCFs ship with the
  ## repository. See data/README.md.

  ## 0. Cohorts. Everything downstream reads these.
  preprocess = c("Preprocess_ICGC_BreastAdenoCA.R",
                 "Preprocess_Breast80.R",
                 "Load_PCAWG_clinical.R"),

  ## 1. Simulations (Figures S1-S3, Tables S1-S4). Days with MCMC on.
  simulations = c("Simulation_main.R",
                  "Simulation_misspec.R"),

  ## 2. The two applications at 2 kb, plus the sensitivity scenarios. Hours.
  applications = c("Application_denovo.R",
                   "Application_refit.R",
                   "Application_denovo_sensitivity.R"),

  ## 3. Supplementary analyses. MAP only.
  extra = c("Application_replicability_80Breast.R",
            "Application_stability_of_covariates.R",
            "Comparison_TensorSignatures.R",
            "HighCN_arm_composition.R"),

  ## 4. Every figure in the paper.
  figures = c("Figure1_mutation_landscape.R",
              "Figure2_goodness_of_fit.R",
              "Reproduce_figures_Application_denovo.R",
              "Reproduce_figures_Application_refit.R",
              "FigureS4_signature_comparison.R",
              "Reproduce_figures_Simulation_misspec.R"),

  ## 5. The LaTeX tables.
  tables = c("Tables_paper.R")
)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) && args[1] %in% c("--list", "-l")) {
  for (nm in names(STAGES)) {
    cat("\n", nm, "\n", sep = "")
    cat(paste0("    ", STAGES[[nm]], collapse = "\n"), "\n")
  }
  quit(save = "no")
}

wanted <- if (length(args)) args else names(STAGES)
unknown <- setdiff(wanted, names(STAGES))
if (length(unknown)) {
  stop("unknown stage(s): ", paste(unknown, collapse = ", "),
       "\n  available: ", paste(names(STAGES), collapse = ", "), call. = FALSE)
}

message("running stages: ", paste(wanted, collapse = ", "))
for (nm in wanted) for (script in STAGES[[nm]]) run(script)
message("\nALL DONE")
