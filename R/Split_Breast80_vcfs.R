################################################################################
# Produces: no figure. Builds data/data_for_breast80/SNP80Breast/
#
# Split the combined CaVEMan table into one VCF per tumour
#
# The Sanger release distributes the 80-cohort substitutions as a SINGLE table -
# Caveman_80sample_Yclean_1Dec15.txt, one row per substitution across all 80
# tumours - but Preprocess_Breast80.R reads per-sample VCFs. This rebuilds
# data/SNP80Breast/ from that download, so the cohort can be reconstructed from
# the public source rather than from the copies committed here.
#
#   wget https://ftp.sanger.ac.uk/pub/cancer/Nik-ZainalEtAl-560BreastGenomes/Caveman_80sample_Yclean_1Dec15.txt
#
# Usage:  Rscript R/Split_Breast80_vcfs.R [table] [outdir]
#
#   table   the combined download. Default data/Caveman_80sample_Yclean_1Dec15.txt
#   outdir  where the VCFs go. Default data/SNP80Breast, which this REFUSES to
#           overwrite - those are the published inputs. Delete them first, or
#           point somewhere else and diff.
#
# The table carries more than the VCFs keep: the gene/transcript annotation
# columns (redundant with the VD tag), the dbSNP columns DS and SNP, and the
# fourteen FLG-* post-hoc filter flags are all dropped, exactly as in the VCFs
# committed here - none of them is declared in those files' headers.
#
# Numbers are reformatted, not copied. The table writes probabilities in
# scientific notation throughout (9.9e-01, 0.0e+00); the VCFs hold what R's
# as.character() makes of them (0.99, 0), which switches to fixed notation only
# when that is no wider. Reproducing the VCFs byte for byte means going through
# as.numeric() rather than passing the strings along.
################################################################################

suppressPackageStartupMessages(library(data.table))

## Run from the repository root, or from R/.
source(if (file.exists("config.R")) "config.R" else "../config.R")

args   <- commandArgs(trailingOnly = TRUE)
TABLE  <- if (length(args) >= 1) args[1] else PATH_BREAST80_CAVEMAN
OUTDIR <- if (length(args) >= 2) args[2] else PATH_BREAST80_SNV

if (!file.exists(TABLE)) {
  stop("no such file: ", TABLE,
       "\n\nDownload it with\n  wget https://ftp.sanger.ac.uk/pub/cancer/",
       "Nik-ZainalEtAl-560BreastGenomes/Caveman_80sample_Yclean_1Dec15.txt",
       "\nSee data/README.md.", call. = FALSE)
}

## These are published inputs; MANIFEST.tsv records their checksums. Rebuilding
## over the top of them would quietly invalidate that.
existing <- list.files(OUTDIR, pattern = "\\.vcf$")
if (length(existing)) {
  stop(length(existing), " VCFs already exist in ", OUTDIR,
       ".\n  Delete them to rebuild, or pass a different output directory:",
       "\n    Rscript R/Split_Breast80_vcfs.R ", TABLE, " /tmp/vcf_check",
       call. = FALSE)
}
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

## ------------------------------------------------------------- the VCF header
## Verbatim from the files this reproduces. Only the #CHROM line varies, and only
## in the two sample names it ends with - normal first, then tumour.
VCF_HEADER <- c(
  "##fileformat=VCFv4.2",
  "##source=CaVEMan_NikZainal_560BreastGenomes",
  "##reference=GRCh37",
  "##FILTER=<ID=PASS,Description=\"All filters passed\">",
  "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">",
  "##INFO=<ID=MP,Number=1,Type=Float,Description=\"Sum of somatic genotype probabilities\">",
  "##INFO=<ID=GP,Number=1,Type=Float,Description=\"Sum of germline genotype probabilities\">",
  "##INFO=<ID=TG,Number=1,Type=String,Description=\"Most probable genotype (CaVEMan)\">",
  "##INFO=<ID=TP,Number=1,Type=Float,Description=\"Probability of most probable genotype\">",
  "##INFO=<ID=SG,Number=1,Type=String,Description=\"2nd most probable genotype\">",
  "##INFO=<ID=SP,Number=1,Type=Float,Description=\"Probability of 2nd most probable genotype\">",
  "##INFO=<ID=VD,Number=1,Type=String,Description=\"Vagrent default annotation\">",
  "##INFO=<ID=VW,Number=1,Type=String,Description=\"Vagrent most deleterious annotation\">",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  "##FORMAT=<ID=FAZ,Number=1,Type=Integer,Description=\"+ strand A reads\">",
  "##FORMAT=<ID=FCZ,Number=1,Type=Integer,Description=\"+ strand C reads\">",
  "##FORMAT=<ID=FGZ,Number=1,Type=Integer,Description=\"+ strand G reads\">",
  "##FORMAT=<ID=FTZ,Number=1,Type=Integer,Description=\"+ strand T reads\">",
  "##FORMAT=<ID=RAZ,Number=1,Type=Integer,Description=\"- strand A reads\">",
  "##FORMAT=<ID=RCZ,Number=1,Type=Integer,Description=\"- strand C reads\">",
  "##FORMAT=<ID=RGZ,Number=1,Type=Integer,Description=\"- strand G reads\">",
  "##FORMAT=<ID=RTZ,Number=1,Type=Integer,Description=\"- strand T reads\">",
  "##FORMAT=<ID=PM,Number=1,Type=Float,Description=\"Proportion of mutant allele\">")

FORMAT_FIELD <- "GT:FAZ:FCZ:FGZ:FTZ:RAZ:RCZ:RGZ:RTZ:PM"
COUNTS <- c("FAZ", "FCZ", "FGZ", "FTZ", "RAZ", "RCZ", "RGZ", "RTZ")

## Scientific notation in, R's default numeric formatting out. Anything that is
## not a number is passed through untouched rather than turned into "NA".
renumber <- function(x) {
  v   <- suppressWarnings(as.numeric(x))
  out <- as.character(v)
  out[is.na(v)] <- x[is.na(v)]
  out
}

## ------------------------------------------------------------------ read once
message("reading ", basename(TABLE))
## Everything as character: the genotype and annotation fields must survive
## verbatim, and the float fields are reformatted explicitly below.
##
## The file opens with a ~64-line two-column legend of ## lines, which fread
## would otherwise try to reconcile with the 62-column table below it. Locate the
## header - the first line on a single # - and skip past the legend by count.
## (skip="#AnalysisProc" does NOT work: it matches the legend's ##AnalysisProc
## line first, on a substring.)
con <- file(TABLE, "r")
preamble <- readLines(con, n = 500L, warn = FALSE)
close(con)
hdr <- which(startsWith(preamble, "#") & !startsWith(preamble, "##"))[1]
if (is.na(hdr)) {
  stop("found no column header in the first 500 lines of ", basename(TABLE),
       "\n  Expected a line starting '#AnalysisProc'.", call. = FALSE)
}
dt <- fread(TABLE, sep = "\t", header = TRUE, colClasses = "character",
            skip = hdr - 1L, showProgress = FALSE)
setnames(dt, sub("^#", "", names(dt)))
message(format(nrow(dt), big.mark = ","), " substitutions | ",
        uniqueN(dt$Sample), " samples")

needed <- c("Sample", "Normal", "VariantID", "Chrom", "Pos", "Ref", "Alt",
            "Qual", "Filter", "DP", "MP", "GP", "TG", "TP", "SG", "SP",
            "VD", "VW",
            paste0(c("GT", COUNTS, "PM"), "-Norm"),
            paste0(c("GT", COUNTS, "PM"), "-Tum"))
missing <- setdiff(needed, names(dt))
if (length(missing)) {
  stop("the table is missing expected column(s): ", paste(missing, collapse = ", "),
       "\n  Is this the 80-sample CaVEMan file?", call. = FALSE)
}

## ------------------------------------------------------------ build the fields
for (j in c("MP", "GP", "TP", "SP", "PM-Norm", "PM-Tum")) {
  set(dt, j = j, value = renumber(dt[[j]]))
}

## The table marks an unannotated variant with "-", VCF with ".". Every record
## Vagrent had nothing to say about differs in exactly this way.
for (j in c("VD", "VW")) set(dt, i = which(dt[[j]] == "-"), j = j, value = ".")

dt[, INFO := paste0("DP=", DP, ";MP=", MP, ";GP=", GP, ";TG=", TG,
                    ";TP=", TP, ";SG=", SG, ";SP=", SP,
                    ";VD=", VD, ";VW=", VW)]

paste_sample <- function(suffix) {
  cols <- paste0(c("GT", COUNTS, "PM"), suffix)
  do.call(paste, c(lapply(cols, function(k) dt[[k]]), sep = ":"))
}
dt[, NORMAL_FIELD := paste_sample("-Norm")]
dt[, TUMOUR_FIELD := paste_sample("-Tum")]

## Lexical on chromosome, numeric on position - "1", "10", "11", ..., "2", ...,
## "X", which is the order the shipped VCFs are in. setorder sorts characters in
## the C locale regardless of the session's, so this does not drift.
dt[, POS_NUM := as.integer(Pos)]
setorder(dt, Chrom, POS_NUM)

dt[, LINE := paste(Chrom, Pos, VariantID, Ref, Alt, Qual, Filter, INFO,
                   FORMAT_FIELD, NORMAL_FIELD, TUMOUR_FIELD, sep = "\t")]

## ------------------------------------------------------------------ write out
samples <- sort(unique(dt$Sample))
message("writing ", length(samples), " VCFs to ", OUTDIR)
for (s in samples) {
  rows   <- dt[Sample == s]
  normal <- unique(rows$Normal)
  if (length(normal) != 1L) {
    stop("sample ", s, " has ", length(normal), " matched normals: ",
         paste(normal, collapse = ", "), call. = FALSE)
  }
  chrom_line <- paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                        "INFO", "FORMAT", normal, s), collapse = "\t")
  writeLines(c(VCF_HEADER, chrom_line, rows$LINE),
             file.path(OUTDIR, paste0(s, ".caveman.vcf")))
}

message(sprintf("done: %d files, %s records",
                length(samples), format(nrow(dt), big.mark = ",")))
