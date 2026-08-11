# Inputs

Not tracked: the ICGC cohort is access-controlled and cannot be redistributed.

Expected here (or wherever `SIGNATUREPPF_DATA` points):

- `Breast80_data.rds.gzip` — Davies et al. (2017), 80 breast WGS, 10 kb bins
- `ICGC_BreastAdenoCA_avg10kb_Mutations_Covariates_Copies.rds.gzip` — ICGC Breast-AdenoCa, 10 kb bins
- `data_for_application/hg19-blacklist.v2.bed`
- `data_for_application/gaps_hg19.bed`
- `data_for_application/Breast-AdenoCa_snp.rds.gzip`
- `data_for_application/20170119_final_consensus_copynumber_donor`

The ChromHMM segmentation (Roadmap E028, breast epithelium) is a public download:

    wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E028_15_coreMarks_dense.bed.gz

Point `SIGNATUREPPF_CHROMHMM` at the unzipped file.
