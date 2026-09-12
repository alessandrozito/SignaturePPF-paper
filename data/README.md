# Inputs

| Folder | Size | Contents |
|---|---|---|
| `data_for_application/` | 336 MB | ICGC mutations, copy number, ChromHMM, masks, and the 18 covariate bigWigs |
| `data_for_breast80/` | 98 MB | 80 CaVEMan VCFs and 80 ASCAT segment tables |
| `data_for_figure1/` | 17 MB | three files from `data_for_application/`, enough to redraw Figure 1 alone |
| `pcawg_clinical/` | 2 MB | PCAWG donor annotation, fetched by `R/Load_PCAWG_clinical.R` |
| `preprocessed/` | — | binned cohorts, written by the `Preprocess_*` scripts. Not tracked. |

Check a copy against what the published results used:

```sh
cd data && md5sum -c <(awk -F'\t' 'NR>1 {print $3"  "$1}' MANIFEST.tsv)
```

To read the data from elsewhere, `export SIGNATUREPPF_DATA=/path/to/data`.

## Sources

### ICGC Breast-AdenoCa — UCSC Xena, PCAWG hub

Open access: the `.nonUS` release covers the non-US ICGC projects, and the 113
donors used here are all BRCA-EU or BRCA-UK.

```sh
wget https://pcawg-hub.s3.us-east-1.amazonaws.com/download/20170119_final_consensus_copynumber_donor
wget https://pcawg-hub.s3.us-east-1.amazonaws.com/download/October_2016_whitelist_2583.snv_mnv_indel.maf.xena.nonUS
```

Browser pages:
[copy number](https://xenabrowser.net/datapages/?dataset=20170119_final_consensus_copynumber_donor&host=https%3A%2F%2Fpcawg.xenahubs.net),
[mutations](https://xenabrowser.net/datapages/?dataset=October_2016_whitelist_2583.snv_mnv_indel.maf.xena.nonUS&host=https%3A%2F%2Fpcawg.xenahubs.net).
`https://pcawg.xenahubs.net/download/<dataset>` serves the same bytes.

The copy-number file is used as downloaded and matches `MANIFEST.tsv` exactly.

`Breast-AdenoCa_snp.rds.gzip` is the Breast-AdenoCa slice of the MAF as a
`GRanges` — 713,855 substitutions, 113 donors, with `tumor`, `sample` and
`channel`. Built by subsetting to those donors, keeping single-base
substitutions, renaming `1` to `chr1`, and assigning each mutation its
trinucleotide channel from hg19. **No script here does this yet**; the object
comes from the predecessor project.

### 80 breast cancers — Nik-Zainal et al. (2016), Sanger

```sh
wget -r -np -nH --cut-dirs=5 -R "index.html*" -P copyNumber80Breast \
  https://ftp.sanger.ac.uk/pub/cancer/Nik-ZainalEtAl-560BreastGenomes/CopyNumber/Breast80/
wget https://ftp.sanger.ac.uk/pub/cancer/Nik-ZainalEtAl-560BreastGenomes/Caveman_80sample_Yclean_1Dec15.txt
```

Copy number is ready to use. The substitutions arrive as one combined table;
`Rscript R/Split_Breast80_vcfs.R` turns it into the 80 per-sample VCFs, byte for
byte. Pass an output directory to compare instead of overwriting:

```sh
Rscript R/Split_Breast80_vcfs.R data/Caveman_80sample_Yclean_1Dec15.txt /tmp/vcf_check
```

### Masks and annotation

```sh
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E028_15_coreMarks_dense.bed.gz
wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg19-blacklist.v2.bed.gz
```

Assembly gaps: UCSC Table Browser, hg19, Mapping and Sequencing → Gap, as BED.

PCAWG clinical: `Rscript R/Load_PCAWG_clinical.R` downloads all four files from
`https://object.genomeinformatics.org/icgc25k-open` and writes the joined table
to `output/PCAWG_clinical/`. Read it with `load_clinical()`. Note that
`donor_survival_time` is empty for every donor in this cohort.

### Covariate tracks

| File | Covariate |
|---|---|
| `gc_content_1kb.bigWig` | GC |
| `Breast-Cancer_Methylation.bigWig` | Methyl |
| `Breast-Cancer_{tissue,cell}_<MARK>_2kb.bigWig` | CTCF, H3K9me3, H3K36me3, H3K27me3, H3K27ac, H3K4me1, H3K4me3 — the two sources are averaged |
| `wgEncodeUwRepliSeqMcf7WaveSignalRep1.bigWig` | RepliTime |
| `GSM920557_hg19_wgEncodeSydhNsomeK562Sig_1kb.bigWig` | NuclOccup |

## preprocessed/

| File | Built by |
|---|---|
| `ICGC_BreastAdenoCA_avg2kb_...rds.gzip` | `Rscript R/Preprocess_ICGC_BreastAdenoCA.R` |
| `ICGC_BreastAdenoCA_avg10kb_...rds.gzip` | `Rscript R/Preprocess_ICGC_BreastAdenoCA.R 10000` |
| `Breast80_data.rds.gzip` | `Rscript R/Preprocess_Breast80.R` |

Each is a list of `gr_Mutations`, `SignalTrack` and `CopyTrack`. The 10 kb file
currently on disk predates the `merge_with_tumor()` fix described in `NOTES.md`,
so rebuilding it will not match byte for byte.
