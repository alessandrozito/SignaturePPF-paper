################################################################################
# Produces: no figure. Builds output/PCAWG_clinical/
#
# Donor-level clinical annotation for the ICGC Breast-AdenoCa cohort
#
# Usage:  Rscript R/Load_PCAWG_clinical.R
#
# -----------------------------------
# The cohort is drawn from two ICGC projects with different recruitment
# criteria, so `project_code` is a coarse receptor-status label that needs no
# extra data:
#
#   BRCA-EU   Breast ER+ HER2- Cancer, EU/UK
#   BRCA-UK   Breast Triple Negative / Lobular Cancer, UK
#
################################################################################

suppressPackageStartupMessages({
  library(readxl)
})

## Run from the repository root, or from R/.
source(if (file.exists("config.R")) "config.R" else "../config.R")
load_functions()

################################################################################
# 1. Fetch, cached
################################################################################
message("== 1. inputs ==")
for (dest in names(PCAWG_REMOTE)) {
  if (file.exists(dest)) {
    message("cached: ", basename(dest))
    next
  }
  url <- paste(PCAWG_BUCKET, PCAWG_REMOTE[[dest]], sep = "/")
  message("downloading: ", basename(dest))
  # Written to a temporary name and moved into place, so an interrupted download
  # leaves nothing that a later run would mistake for a complete cache entry.
  tmp <- paste0(dest, ".part")
  ok <- utils::download.file(url, tmp, mode = "wb", quiet = TRUE)
  if (ok != 0 || !file.exists(tmp) || file.size(tmp) == 0) {
    unlink(tmp)
    stop("could not download ", url, call. = FALSE)
  }
  file.rename(tmp, dest)
}

################################################################################
# 2. Join onto the cohort
################################################################################
message("\n== 2. join ==")
data <- readRDS(PATH_ICGC2KB)
donors <- colnames(data$CopyTrack)

clinical  <- as.data.frame(readxl::read_excel(PATH_PCAWG_CLINICAL))
histology <- as.data.frame(readxl::read_excel(PATH_PCAWG_HISTOLOGY))
sheet     <- utils::read.delim(PATH_PCAWG_SAMPLESHEET, stringsAsFactors = FALSE)

missing <- setdiff(donors, clinical$icgc_donor_id)
if (length(missing)) {
  stop(length(missing), " donor(s) are not in the PCAWG clinical table: ",
       paste(utils::head(missing, 5), collapse = ", "), call. = FALSE)
}

# The Sanger identifier, recovered from the TUMOUR specimen row. A donor has
# several rows here - blood normal, tissue normal, tumour - and only the tumour
# one carries the PDxxxxa label the Sanger releases and the Nik-Zainal
# supplementary tables are keyed by.
tum <- sheet[sheet$dcc_specimen_type == "Primary tumour - solid tissue", ]
pd <- sub("[a-z]+[0-9]*$", "", tum$submitter_sample_id[match(donors, tum$icgc_donor_id)])

ci <- match(donors, clinical$icgc_donor_id)
hi <- match(donors, histology$icgc_donor_id)

num <- function(x) suppressWarnings(as.numeric(x))
out <- data.frame(
  donor        = donors,
  PD           = pd,
  project      = clinical$project_code[ci],
  sex          = clinical$donor_sex[ci],
  age          = num(clinical$donor_age_at_diagnosis[ci]),
  vital        = clinical$donor_vital_status[ci],
  survival     = num(clinical$donor_survival_time[ci]),
  last_followup = num(clinical$donor_interval_of_last_followup[ci]),
  therapy      = clinical$first_therapy_type[ci],
  grade        = histology$tumour_grade[hi],
  stage        = histology$tumour_stage[hi],
  histology    = histology$tumour_histological_type[hi],
  cellularity  = num(histology$percentage_cellularity[hi]),
  stringsAsFactors = FALSE)

# T stage on its own, the part with enough donors behind it to be worth using.
out$T_stage <- ifelse(is.na(out$stage), NA_character_,
                      toupper(sub("^(T[0-9Xx]).*", "\\1", out$stage)))

utils::write.csv(out, PATH_CLINICAL, row.names = FALSE)
saveRDS(out, sub("\\.csv$", ".rds.gzip", PATH_CLINICAL), compress = "gzip")

################################################################################
# 3. Checks
################################################################################
message("\n== 3. checks ==")

# Documented above: if this ever stops being true, the header is out of date.
if (any(!is.na(out$survival))) {
  warning("`donor_survival_time` is now populated for ",
          sum(!is.na(out$survival)), " donor(s). The note at the top of this ",
          "file says it is empty for all of them - revisit it.", call. = FALSE)
}

# The replication analysis compares this cohort against the 80-cancer one, which
# is drawn from the same Sanger releases. Asserting that no donor appears in both
# is what makes "replication" mean replication and not re-fitting.
if (file.exists(PATH_BREAST80)) {
  pd80 <- sub("[a-z]$", "", colnames(readRDS(PATH_BREAST80)$CopyTrack))
  shared <- intersect(out$PD, pd80)
  if (length(shared)) {
    warning("the two breast cohorts share ", length(shared), " donor(s): ",
            paste(utils::head(shared, 5), collapse = ", "),
            ". The replication analysis assumes they are disjoint.",
            call. = FALSE)
  } else {
    message("the ICGC and 80-cancer cohorts are disjoint (0 shared donors).")
  }
}

################################################################################
# 4. Completeness report
################################################################################
message("\n== 4. what is usable ==")
completeness <- data.frame(
  field = setdiff(names(out), "donor"),
  n_present = vapply(setdiff(names(out), "donor"),
                     function(v) sum(!is.na(out[[v]]) & out[[v]] != ""),
                     integer(1)),
  row.names = NULL)
completeness$pct <- round(100 * completeness$n_present / nrow(out))
print(completeness[order(-completeness$n_present), ], row.names = FALSE)

message("\nproject (the receptor-status proxy):")
print(table(out$project))
message("grade:")
print(table(out$grade, useNA = "ifany"))

message("\nwritten: ", PATH_CLINICAL)
