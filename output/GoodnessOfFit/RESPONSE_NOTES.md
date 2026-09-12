# Goodness of fit and patient-level deviations — numbers for the response

All from `Rscript R/GoodnessOfFit.R` on the 2 kb ICGC Breast-AdenoCa de novo MAP
fit (1,391,345 bins, 113 patients, 707,104 mutations). Observed and expected
totals agree to +0.33%.

No posterior predictive check, no patient excluded, nothing refitted.

---

## 0. Why the analysis is built out of effect sizes

The referee's questions are all of the form "is there a systematic residual in
X". At 707k mutations every consistent test answers *yes* for deviations far
below the size that would change a conclusion — the chromosome table below is
the cleanest illustration: 16 of 23 chromosomes are "significant" at BH 0.05
while **no chromosome's count is off by more than 19%, and only one by more than
9%**. So the reported quantities are effect sizes, and the tests that remain are
run at power held equal across slices.

The effect size is the **excess rate error**. If the model's regional rate
carries a multiplicative error with coefficient of variation `cv` on top of
Poisson noise, the dispersion ratio satisfies

```
E[X²] = Σ (e_r + cv²e_r²)/e_r = R + cv²Σe_r    ⟹   D = 1 + cv²·mean(e)
```

so `cv = sqrt((D−1)/mean(e))` is the fractional error in the predicted rate.
`D` is judged against its **exact** Poisson sampling band — for `Y = (X−e)²/e`,
`E[Y] = 1` and `Var(Y) = 2 + 1/e` — not chi-square asymptotics, so it stays valid
at 2 kb where the expected counts are ≪ 1. No minimum expected count is imposed.

---

## 1. The 2 kb count distribution (the referee's own example)

Closed form, not simulated: the counts in disjoint cells are independent Poisson
with known mean, so the expected number of (bin, patient) cells holding exactly
*m* mutations is `Σ Pois(m; λ_bj)`, with exact variance `Σ p_m(1−p_m)`.

| mutations in a 2 kb bin | observed cells | expected cells | obs/exp | share of all mutations |
|---|---|---|---|---|
| 0 | 156,537,609 | 156,523,000 | **1.0001** | — |
| 1 | 668,063 | 694,158 | **0.962** | 94.5% |
| 2 | 13,353 | 5,151 | 2.59 | 3.8% |
| 3 | 1,741 | 92.6 | 18.8 | 0.74% |
| 4 | 575 | 12.0 | 48.0 | 0.33% |
| ≥5 | 644 | 4.96 | 129.8 | 0.68% |

**The statement to make:** the model reproduces the 2 kb count distribution to
within 0.01% for empty bins and 4% for singly-mutated bins, which together carry
**98.3% of the cohort**. It under-predicts bins carrying three or more mutations
by one to two orders of magnitude, and those bins hold **1.7%** of the
mutations. That excess is local clustering — kataegis and hotspots — which a
smooth log-linear intensity does not claim to express, and its location is
identified in §3.

---

## 2. Dispersion across scales, against the covariate-free null

The nested null is the same de novo model with β ≡ 0 (`MAPSolution_noCovariates`,
already fitted for Figure 3): same likelihood, same prior, same optimizer.

| scale | D | rate error | D (no covariates) | rate error (no cov.) | variance explained |
|---|---|---|---|---|---|
| 2 kb | 1.17 | 57.4% | 1.18 | 59.2% | 5.8% |
| 10 kb | 1.37 | 38.4% | 1.45 | 42.2% | 17.2% |
| 100 kb | 2.13 | 21.2% | 3.05 | 28.6% | 45.1% |
| 1 Mb | 6.12 | **14.4%** | 13.90 | 22.9% | **60.3%** |
| 10 Mb | 20.96 | **9.3%** | 49.20 | 14.5% | **58.6%** |

**The statement to make:** the predicted regional rate is accurate to 9% at
10 Mb, 14% at 1 Mb and 21% at 100 kb, and the covariates account for 45–60% of
the excess regional variance at those scales. At 2 kb they account for 6% —
fine-scale variation is dominated by clustering the covariates do not describe,
which is the honest limit of the model and is exactly what §1 measures.

Note that `D` rises with scale while the rate error *falls*. That is not a
contradiction: `D = 1 + cv²·mean(e)` and `mean(e)` grows with the region, so `D`
alone is not comparable across scales. This is the whole reason the rate error
is the reported quantity.

Per mutation class the ordering is the same at every scale; C>T is the worst and
T>G / T>C the best (`dispersion_by_class_and_scale.csv`).

---

## 3. Where the misfit is — regions, chromosomes, mutation types

**Regions (1 Mb).** 707 of 2,867 flagged (BH 0.05), but the misfit is
concentrated: the worst 1% of regions (29) carry **19.6%** of it, and removing
the worst 5% takes `D` from 6.12 to 3.77.

The largest departures are interpretable, and they run in both directions:

| region | observed | expected | resid | what it is |
|---|---|---|---|---|
| chr14:107.0–107.3 Mb | 89 | 11.1 | +23.4 | IGH locus (AID / kataegis) |
| chr6:126–127 Mb | 574 | 230.4 | +22.6 | 100% callable — unexplained |
| chr14:106–107 Mb | 250 | 88.9 | +17.1 | IGH locus |
| chr9:39–40 Mb | 7 | 174.5 | **−12.7** | 99.7% callable — unexplained |
| chr9:40.0–41 Mb | 5 | 137.5 | **−11.3** | 99.2% callable — unexplained |
| chrX:91 Mb | 71 | 215.6 | −9.9 | 100% callable — unexplained |

**These are NOT masking artefacts, and it was checked.** The preprocessing
subtracts assembly gaps and the ENCODE blacklist from each bin's usable width,
and **zero of the 707,104 mutations** fall inside either mask; the fifteen worst
megabases are 96.5–100% usable sequence against a cohort median of 100%. So the
large deficits sit in fully callable sequence: the model predicts 175 mutations
at chr9:39 Mb and 7 are observed. That is a real departure, not under-calling.

Regions are labelled by coordinate, with the immunoglobulin loci named. **No gene
symbols**: the longest gene overlapping a megabase is a coordinate lookup, and
printed on a residual plot it reads as a claim that the model fails at that gene,
which nothing supports.

**Chromosomes.** The one cut where the total is not fitted away — the model has
no chromosome term. Fold change spans **0.915 to 1.185 across all 23**; only chrX
(+18.5%) exceeds ±10%, and it is the expected one. 16 of 23 reach BH 0.05, which
is the large-*n* point again.

**Mutation types.** Flagged 1 Mb regions, and the share of that class's misfit
they carry: C>T 311 regions / 60.3%, C>A 156 / 36.7%, T>A 51 / 23.9%, C>G 35 /
15.8%, T>C 8 / 8.8%, T>G 1 / 0.8%. The localised misfit is a C>T phenomenon,
consistent with the IGH/kataegis regions above; T>C and T>G show essentially
none.

---

## 4. Patient-level deviations — the outlier question

**Two tests, because two questions.** The time-rescaling test on ALL of a
patient's mutations (`p_full`) is what says whether that patient is
misspecified. The same test after binomial thinning to a common event count
(`p_thin`) is *only* for asking whether misfit tracks burden, where the
burden-driven power difference is the confounder.

```
rejected on all their mutations : 86 of 113
rejected after thinning to 1,190 : 34 of 113
```

The gap is not a contradiction. Thinning does not shrink a deviation, it raises
the bar: a spike of *k* events moves the empirical CDF by about `k/n`, which
thinning leaves alone, while the critical distance grows like `1/sqrt(n)`. For
DO1020, the KS distance stays 0.0552 while `D*` moves 0.0134 → 0.0394. Correct
for a power-matched comparison; wrong for detection. **Report 86/113.**

**Association with burden** (thinned arm, and the burden-invariant effect sizes):

| quantity | burden-invariant | Spearman ρ | p |
|---|---|---|---|
| excess rate error (cv) | yes | **−0.200** | 0.034 |
| −log₁₀ equal-power p | yes | **−0.243** | 0.009 |
| share of mutations in flagged regions | yes | **−0.188** | 0.046 |
| χ² in worst 0.1% of regions | yes | +0.036 | 0.71 |
| dispersion ratio D | **no** | +0.470 | 1.5e-07 |
| KS distance at full n | **no** | −0.396 | 1.4e-05 |

The last two rows are mechanically confounded and are shown with the reason
attached: `D = 1 + cv²·mean(e)` with `mean(e) ∝ burden`, and the KS distance at
full `n` grows like √n under any fixed misspecification. Every burden-invariant
quantity gives zero or a mild **negative** association — hypermutated patients
are, if anything, fitted slightly better in fractional terms. **There is no case
for excluding them**, which is worth saying explicitly rather than performing the
exclusion and reporting that it changed nothing.

**Which patients, and why.** The volcano panel plots KS distance against
significance, one point per patient, sized by burden, with focal patients
outlined. The extreme is **DO1020**: KS distance 0.0552 against a critical 0.0134,
`padj = 5e-26`, rate error 1.67 versus a median of 0.46, and **87% of its
departure inside 0.1% of its megabases** — one patient with one hotspot
(chr6:126 Mb, 363 observed against 5.5 expected), not a globally misfitted
patient. It is also the sample the covariate-stability analysis excluded by hand,
which the diagnostic recovers without being told.

The (patient, megabase) cell volcano in the supplement adds the other arm:
**DO218404 at chr22:49–50 Mb has 2 observed against 154 expected, and 4 against
97** — a large *deficit* in one patient, which reads as uncalled copy-number loss
rather than a failure of the intensity model.

## Outputs

```
output/GoodnessOfFit/
  count_distribution_2kb.csv          §1, all classes
  dispersion_by_class_and_scale.csv   §2, 7 classes x 5 scales
  covariate_gain.csv                  §2, against the covariate-free null
  regions_1Mb_all.csv                 §3, every 1 Mb region
  regions_1Mb_by_class.csv            §3, per mutation class
  top_misfit_regions.csv              §3, worst 25, with usable-sequence fraction
  chromosomes.csv                     §3, per chromosome
  flagged_regions_by_class.csv        §3, class summary
  misfit_concentration.csv            §3, Lorenz summary
  patient_deviations.csv              §4, one row per patient
  burden_association.csv              §4, the table above
  Table_goodness_of_fit.csv           the summary table

figures/
  GoF_count_distribution_2kb.pdf      §1
  GoF_dispersion_by_scale.pdf         §2, D and rate error, by class and scale
  GoF_residuals_along_genome.pdf      §3, 1 Mb residuals, top loci named
  GoF_misfit_concentration.pdf        §3, Lorenz curve
  Figure3_goodness_of_fit.pdf         the paper figure
  FigureS_residuals_by_mutation_class.pdf   §3, per macro class
  FigureS_patient_region_volcano.pdf  §4, (patient, Mb) cells
  FigureS_patient_qq_equal_power.pdf  §4, the thinned arm
```

## Not done, and available if wanted

A **leave-out refit** dropping the hypermutated patients. §4 argues it is not
needed — the association it would probe is absent once power and `mean(e)` are
divided out — but if the referee insists on seeing it, the cheap version is one
MAP restart warm-started from the current mode with those patients removed, plus
a β overlay against the full-cohort fit. A few hours at 2 kb; minutes at 10 kb.
