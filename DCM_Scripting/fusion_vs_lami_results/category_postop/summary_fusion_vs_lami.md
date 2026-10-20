### Laminectomy only vs laminectomy + fusion (PostOp, no preload)

**H0:** adding fusion does not change % cord volume above MPS = 0.015; the paired difference d = fusion − lami is symmetric about 0. Flexion and Extension are analysed separately (not pooled). Category source: PostOp.

Primary test: exact Wilcoxon signed-rank with Pratt zero handling. Flexion and Extension are separate hypotheses, each tested once, so no multiplicity adjustment is applied. The other thresholds and all other tests are supporting or sensitivity analyses.

#### Headline results

| Outcome | Condition | Median difference, fusion - lami (pp) (95% CI) | Wilcoxon signed-rank p |
|---|---|---|---|
| Cord, MPS ≥ 0.015 | Flexion | -10.03 (-23.75 to +0.00) | 0.004 |
| Cord, MPS ≥ 0.015 | Extension | +0.00 (-0.06 to +0.19) | 1.000 |
| IVD, global T95 | Flexion | +1.32 (+0.52 to +2.10) | 0.034 |
| IVD, global T95 | Extension | +2.12 (+0.44 to +4.07) | 0.027 |

Median difference = median of the 12 paired differences (fusion - lami), in percentage points of volume above the threshold. The 95% CI is the exact distribution-free interval for the median (order statistics, binomial; coverage at least 95%, 96.1% for 12 pairs). p = exact Wilcoxon signed-rank test (Pratt zero handling). Flexion and extension are separate hypotheses, each tested once.

#### Flexion (threshold = 0.015)

**12 pairs; 9 differ, 3 are tied at zero (both arms have no volume above threshold).**

| Pair type | n |
|---|---|
| Both arms zero | 3 |
| Lami only above threshold | 2 |
| Fusion only above threshold | 0 |
| Both above threshold | 7 |

| Test | Estimate | p | Note |
|---|---|---|---|
| **Wilcoxon signed-rank, exact, Pratt (primary)** | W+ = 0.0 | **0.004** | m = 9 non-zero |
| Permutation test, mean difference | -11.17 pp (95% CI -18.48 to -3.96) | 0.004 | exact, 2^9 sign patterns |
| Hodges-Lehmann median difference | -11.38 pp | | |
| Sign test (\|d\| > 0.01 pp) | lower 9, higher 0, tied 3 | 0.004 | |
| KR LME: surgery + (1 \| patient) | -11.17 pp (95% CI -18.37 to -3.97) | 0.006 | Shapiro-Wilk (LME residuals) W = 0.864, p = 0.004 (deviates from normality) |
| paired t-test (LME singular): surgery + category + (1 \| patient) | -11.17 pp (95% CI -18.37 to -3.97) | 0.006 | Shapiro-Wilk (paired differences) W = 0.886, p = 0.103 |
| Exploratory: surgery × category (KR F-test) | | 0.178 | singular fit |

Secondary, **confounded** (surgeons chose who got fusion): factual-only fusion (n=4) vs lami (n=8), mean difference -12.75 pp, exact permutation p = 0.107 (mean), 0.137 (ranks).

#### Extension (threshold = 0.015)

**12 pairs; 7 differ, 5 are tied at zero (both arms have no volume above threshold).**

| Pair type | n |
|---|---|
| Both arms zero | 5 |
| Lami only above threshold | 1 |
| Fusion only above threshold | 1 |
| Both above threshold | 5 |

| Test | Estimate | p | Note |
|---|---|---|---|
| **Wilcoxon signed-rank, exact, Pratt (primary)** | W+ = 31.0 | **1.000** | m = 7 non-zero |
| Permutation test, mean difference | +1.51 pp (95% CI -0.09 to 4.55) | 0.438 | exact, 2^7 sign patterns |
| Hodges-Lehmann median difference | +0.00 pp | | |
| Sign test (\|d\| > 0.01 pp) | lower 3, higher 3, tied 6 | 1.000 | |
| KR LME: surgery + (1 \| patient) | +1.51 pp (95% CI -1.40 to 4.41) | 0.278 | Shapiro-Wilk (LME residuals) W = 0.478, p = <0.001 (deviates from normality) |
| KR LME: surgery + category + (1 \| patient) | +1.51 pp (95% CI -1.40 to 4.41) | 0.278 | Shapiro-Wilk (LME residuals) W = 0.581, p = <0.001 (deviates from normality) |
| Exploratory: surgery × category (KR F-test) | | 0.626 |  |

Secondary, **confounded** (surgeons chose who got fusion): factual-only fusion (n=4) vs lami (n=8), mean difference +0.23 pp, exact permutation p = 0.675 (mean), 0.582 (ranks).

#### Sensitivity to threshold

The thresholds are one question measured four ways, not four tests; 0.015 was fixed as primary in advance. Higher thresholds produce more tied-at-zero pairs.

| Condition | Threshold | Differ / pairs | Mean diff (pp) | Signed-rank p | Permutation p | Sign test p |
|---|---|---|---|---|---|---|
| Flexion | 0.01 | 9/12 | -15.12 | 0.004 | 0.004 | 0.004 |
| Flexion | 0.015 | 9/12 | -11.17 | 0.004 | 0.004 | 0.004 |
| Flexion | 0.02 | 9/12 | -8.59 | 0.004 | 0.004 | 0.004 |
| Flexion | 0.025 | 9/12 | -6.72 | 0.004 | 0.004 | 0.004 |
| Extension | 0.01 | 8/12 | +1.73 | 0.844 | 0.438 | 0.727 |
| Extension | 0.015 | 7/12 | +1.51 | 1.000 | 0.438 | 1.000 |
| Extension | 0.02 | 7/12 | +1.39 | 1.000 | 0.438 | 1.000 |
| Extension | 0.025 | 7/12 | +1.26 | 0.688 | 0.422 | 1.000 |

#### Interpretation notes

- Pairs tied at zero stay in the mean but carry no sign information, so the exact tests can only reach p = 2/2^m (m = number of differing pairs). Subgroups of 3 to 5 patients can never reach p < 0.05; category effects are shown descriptively (paired-difference plot) and the interaction is exploratory.
- These are deterministic simulations: p-values describe this cohort, not new patients.
- Exactly-zero strain in many PostOp models is physically plausible (a decompressed cord with no preload need not touch bone or ligament) but has not been confirmed for every job.
