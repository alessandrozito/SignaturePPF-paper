################################################################################
# Produces: Tables S1, S2, S3, S4 and S5
#
# LaTeX tables for the manuscript
#
# The MAP timing table and the MCMC ESS table for BOTH simulation studies, from
# the same scored columns, plus the de novo prior-sensitivity table. Each is
# written whole, so the .tex files can be \input directly.
#
#   Rscript R/Tables_paper.R
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
})

## Run from the repository root, or from R/.
source(if (file.exists("config.R")) "config.R" else "../config.R")

## ------------------------------------------------------------- formatting
# mean on the first line, sd in \scriptsize underneath, as in the manuscript.
fmt <- function(x, d = 1) {
  if (all(is.na(x))) return("---")
  formatC(mean(x, na.rm = TRUE), format = "f", digits = d, big.mark = "")
}
fmt_sd <- function(x, d = 1) {
  if (all(is.na(x)) || sum(!is.na(x)) < 2) return("---")
  sprintf("{\\scriptsize (%s)}",
          formatC(stats::sd(x, na.rm = TRUE), format = "f", digits = d))
}
# Runtimes span four orders of magnitude between CompNMF and the MCMC.
time_digits <- function(x) if (mean(x, na.rm = TRUE) < 0.1) 3 else 1

# Wrap a body of rows in the manuscript's table environment.
wrap_table <- function(body, caption, label, colspec, header,
                       size = "\\footnotesize") {
  c("\\begin{table}[th]",
    paste0("\\caption{", caption, "}"),
    "\\centering",
    size,
    "\\begin{adjustbox}{max width=1\\textwidth,center}",
    paste0("\\begin{tabular}{", colspec, "}"),
    "\\hline",
    header,
    "\\hline",
    body,
    "\\end{tabular}",
    "\\end{adjustbox}",
    paste0("\\label{", label, "}"),
    "\\end{table}")
}

cell <- function(df, col, d) {
  x <- df[[col]]
  c(fmt(x, d), fmt_sd(x, d))
}

## ============================================================ MAIN STUDY
main <- utils::read.delim(file.path(DIR_SIM_MAIN, "simulation_results.tsv")) |>
  mutate(Scen = ifelse(.data$Scenario == "Scenario_A_indep", "A", "B"))

MAP_ROWS <- c(map_TrueCovs = "\\ref{simmod:add_MAPTrue}) & MAP, true $\\bx$",
              CompNMFBase  = "\\ref{simmod:add_compNMF}) & MAP, CompNMF",
              map_CopyOnly = "\\ref{simmod:add_noCovs}) & MAP, no $\\bx$",
              map_Full     = "\\ref{simmod:add_MAP}) & MAP, all $\\bx$")

map_table_main <- function() {
  out <- character()
  for (m in names(MAP_ROWS)) {
    d <- lapply(c("A", "B"), function(s) main[main$model == m & main$Scen == s, ])
    td <- time_digits(main$time[main$model == m])
    v <- lapply(d, function(x) c(cell(x, "time", td), cell(x, "iter", 1)))
    out <- c(out,
      sprintf("%s\n  & %s  & %s\n  & %s  & %s \\\\",
              MAP_ROWS[[m]], v[[1]][1], v[[1]][3], v[[2]][1], v[[2]][3]),
      sprintf("& & %s & %s\n & %s & %s \\\\",
              v[[1]][2], v[[1]][4], v[[2]][2], v[[2]][4]),
      "\\hline")
  }
  paste(out, collapse = "\n")
}

MCMC_ROWS <- c(mcmc_Full        = c("\\ref{simmod:add_MCMC}) MCMC,", "all $\\bx$, $\\Delta_b=100$"),
               mcmc_TrueCovs200 = c("\\ref{simmod:add_MCMC200}) MCMC,", "all $\\bx$, $\\Delta_b=200$"),
               mcmc_TrueCovs500 = c("\\ref{simmod:add_MCMC500}) MCMC,", "all $\\bx$, $\\Delta_b=500$"))
MCMC_COLS <- c("time", "effectiveSigs", "effectiveTheta", "effectiveBetas",
               "effectiveMu", "effectiveSigma2", "effectiveLogPost")

mcmc_table_main <- function() {
  models <- c("mcmc_Full", "mcmc_TrueCovs200", "mcmc_TrueCovs500")
  labs <- list(c("\\ref{simmod:add_MCMC}) MCMC,", "all $\\bx$, $\\Delta_b=100$"),
               c("\\ref{simmod:add_MCMC200}) MCMC,", "all $\\bx$, $\\Delta_b=200$"),
               c("\\ref{simmod:add_MCMC500}) MCMC,", "all $\\bx$, $\\Delta_b=500$"))
  out <- character()
  for (i in seq_along(models)) {
    m <- models[i]
    dA <- main[main$model == m & main$Scen == "A", ]
    dB <- main[main$model == m & main$Scen == "B", ]
    means <- sds <- character()
    for (cl in MCMC_COLS) {
      d <- if (cl == "time") 0 else 0
      means <- c(means, fmt(dA[[cl]], d), fmt(dB[[cl]], d))
      sds <- c(sds, fmt_sd(dA[[cl]], d), fmt_sd(dB[[cl]], d))
    }
    out <- c(out,
      sprintf("%s\n& %s \\\\", labs[[i]][1], paste(means, collapse = " & ")),
      sprintf("%s\n& %s \\\\", labs[[i]][2], paste(sds, collapse = " & ")),
      "\\hline")
  }
  paste(out, collapse = "\n")
}

## =================================================== MISSPECIFICATION STUDY
mis_file <- file.path(DIR_SIM_MISSPEC, "all_results.csv")
SCEN <- c(S0_baseline = "S0", S1_epigenome_v025 = "S1", S2_epigenome_v1 = "S2",
          S3_epigenome_v4 = "S3", S4_hotspots = "S4", S5_cn_noise = "S5",
          S6_opportunity = "S6")
SCEN_DESC <- c(S0_baseline = "correctly specified",
               S1_epigenome_v025 = "$v^2 = 0.25$",
               S2_epigenome_v1 = "$v^2 = 1$",
               S3_epigenome_v4 = "$v^2 = 4$",
               S4_hotspots = "$+$ hotspots",
               S5_cn_noise = "$+$ CN noise",
               S6_opportunity = "$+$ opportunity")

map_table_misspec <- function(mis) {
  models <- c(CompNMF = "CompressiveNMF", SignatureAnalyzer = "SignatureAnalyzer",
              PPF_map = "PPF, MAP")
  out <- character()
  for (sc in names(SCEN)) {
    means <- sds <- character()
    for (m in names(models)) {
      d <- mis[mis$model == m & mis$Scenario == sc, ]
      td <- time_digits(mis$time[mis$model == m])
      means <- c(means, fmt(d$time, td), fmt(d$iter, 1))
      sds <- c(sds, fmt_sd(d$time, td), fmt_sd(d$iter, 1))
    }
    out <- c(out,
      sprintf("%s & %s & %s \\\\", SCEN[[sc]], SCEN_DESC[[sc]],
              paste(means, collapse = " & ")),
      sprintf(" & & %s \\\\", paste(sds, collapse = " & ")),
      "\\hline")
  }
  paste(out, collapse = "\n")
}

mcmc_table_misspec <- function(mis) {
  out <- character()
  for (sc in names(SCEN)) {
    d <- mis[mis$model == "PPF_mcmc" & mis$Scenario == sc, ]
    means <- sds <- character()
    for (cl in MCMC_COLS) {
      means <- c(means, fmt(d[[cl]], 0))
      sds <- c(sds, fmt_sd(d[[cl]], 0))
    }
    out <- c(out,
      sprintf("%s & %s & %s \\\\", SCEN[[sc]], SCEN_DESC[[sc]],
              paste(means, collapse = " & ")),
      sprintf(" & & %s \\\\", paste(sds, collapse = " & ")),
      "\\hline")
  }
  paste(out, collapse = "\n")
}

## ----------------------------------------------------------------- write
HDR_MAP_MAIN <- c(
  " \\multicolumn{2}{|c|}{\\textsc{Method}} & \\multicolumn{2}{c|}{\\textsc{Scenario A}} & \\multicolumn{2}{c|}{\\textsc{Scenario B}} \\\\",
  "& & \\textsc{Time (min)} & \\textsc{Iterations} & \\textsc{Time (min)} & \\textsc{Iterations} \\\\")

HDR_MCMC_MAIN <- c(
  "\\multicolumn{1}{|c|}{\\textsc{Method}}",
  "& \\multicolumn{2}{c|}{\\textsc{Time (min)}}",
  "& \\multicolumn{2}{c|}{\\textsc{ESS} $R$}",
  "& \\multicolumn{2}{c|}{\\textsc{ESS} $\\Theta$}",
  "& \\multicolumn{2}{c|}{\\textsc{ESS} $B$}",
  "& \\multicolumn{2}{c|}{\\textsc{ESS} $\\mu$}",
  "& \\multicolumn{2}{c|}{\\textsc{ESS} $\\sigma^2$}",
  "& \\multicolumn{2}{c|}{\\textsc{ESS} logPost}\\\\",
  "& \\textsc{a} & \\textsc{b} & \\textsc{a} & \\textsc{b} & \\textsc{a} & \\textsc{b} & \\textsc{a} & \\textsc{b} & \\textsc{a} & \\textsc{b} & \\textsc{a} & \\textsc{b} & \\textsc{a} & \\textsc{b} \\\\")

HDR_MAP_MIS <- c(
  "\\multicolumn{2}{|c|}{\\textsc{Scenario}} & \\multicolumn{2}{c|}{\\textsc{CompressiveNMF}} & \\multicolumn{2}{c|}{\\textsc{SignatureAnalyzer}} & \\multicolumn{2}{c|}{\\textsc{PPF, MAP}} \\\\",
  "& & \\textsc{Time (min)} & \\textsc{Iter.} & \\textsc{Time (min)} & \\textsc{Iter.} & \\textsc{Time (min)} & \\textsc{Iter.} \\\\")

HDR_MCMC_MIS <- paste(
  "\\multicolumn{2}{|c|}{\\textsc{Scenario}} & \\textsc{Time (min)}",
  "& \\textsc{ESS} $R$ & \\textsc{ESS} $\\Theta$ & \\textsc{ESS} $B$",
  "& \\textsc{ESS} $\\mu$ & \\textsc{ESS} $\\sigma^2$ & \\textsc{ESS} logPost \\\\")

CAP_MAP_MAIN <- paste(
  "Computation time and number of iterations for MAP-based methods.",
  "Shown are averages across 20 randomly generated datasets, with standard",
  "deviation in parentheses beneath each entry.")
CAP_MCMC_MAIN <- paste(
  "Computation time and effective sample sizes (ESS) for MCMC. Shown are",
  "averages across 20 datasets, with standard deviations in parentheses.",
  "ESS for $R$, $\\Theta$, $B$, $\\mu$ and $\\sigma^2$ is calculated over the",
  "last 1500 draws; the log-posterior is recorded once every ten iterations,",
  "so its ESS is calculated over the 150 recorded values.")
CAP_MAP_MIS <- paste(
  "Computation time and number of iterations for MAP-based methods under",
  "misspecification. Shown are averages across 20 datasets per scenario, with",
  "standard deviations in parentheses beneath each entry. SignatureAnalyzer",
  "reports no iteration count.")
CAP_MCMC_MIS <- paste(
  "Computation time and effective sample sizes (ESS) of the PPF sampler under",
  "misspecification. Shown are averages across 20 datasets per scenario, with",
  "standard deviations in parentheses. ESS for $R$, $\\Theta$, $B$, $\\mu$ and",
  "$\\sigma^2$ is calculated over the last 1500 draws; the log-posterior is",
  "recorded once every ten iterations, so its ESS is calculated over the 150",
  "recorded values.")

w <- function(lines, path) { writeLines(lines, path); message("wrote ", path) }
show <- function(title, lines) cat("\n===== ", title, " =====\n",
                                   paste(lines, collapse = "\n"), "\n", sep = "")

t1 <- wrap_table(map_table_main(), CAP_MAP_MAIN, "tab:sim_time_results",
                 "|cl|cc|cc|", HDR_MAP_MAIN)
t2 <- wrap_table(mcmc_table_main(), CAP_MCMC_MAIN, "tab:mcmc_results",
                 "|l|cc|cc|cc|cc|cc|cc|cc|", HDR_MCMC_MAIN, size = "\\small")
w(t1, file.path(DIR_SIM_MAIN, "table_MAP_timing.tex"))
w(t2, file.path(DIR_SIM_MAIN, "table_MCMC_ess.tex"))
show("MAIN: MAP timing", t1)
show("MAIN: MCMC ESS", t2)

if (file.exists(mis_file)) {
  mis <- read.csv(mis_file)
  if (!"iter" %in% names(mis)) {
    message("\nall_results.csv has no `iter`/ESS columns yet - re-run: ",
            "Rscript R/Simulation_misspec.R score")
  } else {
    t3 <- wrap_table(map_table_misspec(mis), CAP_MAP_MIS, "tab:sim_time_misspec",
                     "|cl|cc|cc|cc|", HDR_MAP_MIS)
    t4 <- wrap_table(mcmc_table_misspec(mis), CAP_MCMC_MIS, "tab:mcmc_misspec",
                     "|cl|c|cccccc|", HDR_MCMC_MIS, size = "\\small")
    w(t3, file.path(DIR_SIM_MISSPEC, "table_MAP_timing.tex"))
    w(t4, file.path(DIR_SIM_MISSPEC, "table_MCMC_ess.tex"))
    show("MISSPEC: MAP timing", t3)
    show("MISSPEC: MCMC ESS", t4)
  }
}


## ============================================ APPLICATION: PRIOR SENSITIVITY
# One fit per row, not twenty replicates, so the cells are single values and
# there is no standard-deviation line underneath.
sens_file <- file.path(DIR_SENSITIVITY, "sensitivity_summary.csv")

sensitivity_table <- function(sens) {
  lab <- ifelse(sens$scenario == "reference", "Reference",
                sprintf("\\texttt{%s}", gsub("_", "\\\\_", sens$scenario)))
  num <- function(x, d) formatC(x, format = "f", digits = d)
  rows <- sprintf(
    "%s & %d & %d & %s & %d & %d & %s & %s & %d & %s & %s \\\\",
    lab, sens$Kmax, sens$c0, formatC(sens$d0, format = "g"),
    sens$n_selected, sens$matched, num(sens$mean_cosine, 3),
    num(sens$min_mu, 1), sens$iterations, num(sens$sec_per_iter, 1),
    num(sens$hours, 2))
  paste(c(rbind(rows, "\\hline")), collapse = "\n")
}

HDR_SENS <- c(
  "\\multicolumn{1}{|c|}{\\textsc{Scenario}} & \\multicolumn{3}{c|}{\\textsc{Setting}} & \\multicolumn{4}{c|}{\\textsc{Signatures}} & \\multicolumn{3}{c|}{\\textsc{Computation}} \\\\",
  "& $K_{\\max}$ & $c_0$ & $d_0$ & \\textsc{Selected} & \\textsc{Recovered} & \\textsc{Mean cos.} & $\\min\\mu$ & \\textsc{Iterations} & \\textsc{Sec./iter.} & \\textsc{Hours} \\\\")

CAP_SENS <- paste(
  "Sensitivity of the \\emph{de novo} fit to the maximum number of signatures",
  "$K_{\\max}$ and to the prior on $\\mu$. Each row is a single fit.",
  "\\textsc{Selected} is the number of signatures retained after pruning,",
  "\\textsc{Recovered} the number of the ten reference signatures whose",
  "one-to-one partner attains a cosine similarity of at least $0.9$, and",
  "\\textsc{Mean cos.} the average cosine similarity over all matched pairs;",
  "the reference row is matched against itself and so attains ten and one by",
  "construction.")

if (file.exists(sens_file)) {
  sens <- read.csv(sens_file)
  t5 <- wrap_table(sensitivity_table(sens), CAP_SENS, "tab:sensitivity_denovo",
                   "|l|ccc|cccc|ccc|", HDR_SENS)
  w(t5, file.path(DIR_SENSITIVITY, "table_sensitivity.tex"))
  show("APPLICATION: prior sensitivity", t5)
} else {
  message("no ", sens_file, " - run R/Reproduce_figures_Application_denovo.R")
}
