# Inputs

Everything the analyses read lives here, so the pipeline depends on nothing
outside this repository.

**The files themselves are not tracked** — the ICGC cohort is access-controlled
and cannot be redistributed, and the rest are large public downloads.
`MANIFEST.tsv` *is* tracked and records the size and MD5 of every file the
published results were computed from, so a rebuilt copy can be checked against
the originals:

```sh
cd data && md5sum -c <(awk 'NR>1 {print $3"  "$1}' MANIFEST.tsv)
```

## What is here

| File | Source |
|---|---|
| `Breast80_data.rds.gzip` | Davies et al. (2017), 80 breast WGS, preprocessed to 10 kb bins with 11 covariates |
| `ICGC_BreastAdenoCA_avg10kb_Mutations_Covariates_Copies.rds.gzip` | ICGC Breast-AdenoCa, same 10 kb grid and covariates |
| `ICGC_BreastAdenoCA_avg2kb_Mutations_Covariates_Copies.rds.gzip` | the same cohort at 2 kb, what the applications use — **built by `R/Preprocess_ICGC_BreastAdenoCA.R`**, not shipped |
| `Breast-AdenoCa_snp.rds.gzip` | ICGC Breast-AdenoCa SNV calls, unbinned |
| `20170119_final_consensus_copynumber_donor` | PCAWG consensus copy number, per donor |
| `E028_15_coreMarks_dense.bed` | Roadmap ChromHMM 15-state segmentation, breast epithelium |
| `hg19-blacklist.v2.bed` | ENCODE blacklist v2, hg19 |
| `gaps_hg19.bed` | UCSC assembly gaps, hg19 |
| `SNP80Breast/` | 80 CaVEMan VCFs, one per tumour (Davies et al. 2017) |
| `copyNumber80Breast/` | 80 ASCAT segment tables, one per tumour |

`Breast80_data.rds.gzip` is built from the last two by
`Rscript R/Preprocess_Breast80.R`. The two directories are recorded in
`MANIFEST.tsv` with a trailing slash; their checksum is the md5 of the sorted
per-file md5s, so it is stable regardless of directory listing order.

### Covariate tracks

The eleven genomic covariates, as the bigWigs they are binned from:

| File(s) | Covariate |
|---|---|
| `gc_content_1kb.bigWig` | GC |
| `Breast-Cancer_Methylation.bigWig` | Methyl |
| `Breast-Cancer_{tissue,cell}_<MARK>_2kb.bigWig` | CTCF, H3K9me3, H3K36me3, H3K27me3, H3K27ac, H3K4me1, H3K4me3 — tissue and cell line are averaged |
| `wgEncodeUwRepliSeqMcf7WaveSignalRep1.bigWig` | RepliTime |
| `GSM920557_hg19_wgEncodeSydhNsomeK562Sig_1kb.bigWig` | NuclOccup |

The `*_data.rds.gzip` / `*_Covariates_Copies.rds.gzip` objects are the
preprocessed cohorts: lists of `gr_Mutations`, `SignalTrack` and `CopyTrack`. The
2 kb one is built from everything above by
`Rscript R/Preprocess_ICGC_BreastAdenoCA.R`. The 10 kb one shipped here was
built by the predecessor project's loader, before the `merge_with_tumor()` fix
described in the top-level README, so rebuilding it will not reproduce it byte
for byte.

## Obtaining the public files

```sh
# Roadmap ChromHMM, breast epithelium (E028)
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E028_15_coreMarks_dense.bed.gz
gunzip E028_15_coreMarks_dense.bed.gz

# ENCODE blacklist v2
wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg19-blacklist.v2.bed.gz
gunzip hg19-blacklist.v2.bed.gz
```

Assembly gaps come from the UCSC Table Browser (hg19, Mapping and Sequencing →
Gap), exported as BED.

## Controlled-access files

The ICGC/PCAWG mutation calls and consensus copy number require an approved DACO
application through the [ICGC Data
Portal](https://dcc.icgc.org/). They cannot be
redistributed here. The Davies et al. cohort is available from the accompanying
publication's data availability statement.

## Using a different location

If the data has to live elsewhere — a scratch volume on a cluster, say — set

```sh
export SIGNATUREPPF_DATA=/path/to/data
```

and the scripts will read from there instead. Nothing else needs changing.
