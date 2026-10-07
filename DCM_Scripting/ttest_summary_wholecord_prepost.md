# Whole cord - PreOp (no preload) vs PostOp (threshold = 0.02)

Paired t-test, not a mixed model: the random-intercept variance here is genuinely 0 (confirmed via 5 optimizers + profile-likelihood CI), since patients don't keep a consistent PreOp/PostOp rank - unlike the GM/WM models, where it's essential (ICC 0.87-0.97). Flexion/Extension kept as separate tests, not pooled.

## Whole cord (Flexion): PreOp (no preload) vs PostOp

$$H_0:\ \mu_{\Delta} = 0, \quad \Delta_i = \text{pct\_above}_{i,\text{PostOp}} - \text{pct\_above}_{i,\text{PreOp}}$$

No difference in % of whole cord volume above MPS $=0.02$ between PreOp (no preload) and PostOp, within Flexion. Paired t-test on each patient's PostOp-minus-PreOp difference (Flexion only) - see the note above on why this isn't a mixed model.

| Term | Estimate | CI_low | CI_high | df | t | Pr(>\|t\|) |
|---|---|---|---|---|---|---|
| PostOp - PreOp (no preload) | -9.713 | -20.28 | 0.8492 | 11 | -2.024 | 0.06795 |

**Shapiro-Wilk (paired differences)**: W = 0.9483, p = 0.6117

## Whole cord (Extension): PreOp (no preload) vs PostOp

$$H_0:\ \mu_{\Delta} = 0, \quad \Delta_i = \text{pct\_above}_{i,\text{PostOp}} - \text{pct\_above}_{i,\text{PreOp}}$$

No difference in % of whole cord volume above MPS $=0.02$ between PreOp (no preload) and PostOp, within Extension. Paired t-test on each patient's PostOp-minus-PreOp difference (Extension only) - see the note above on why this isn't a mixed model.

| Term | Estimate | CI_low | CI_high | df | t | Pr(>\|t\|) |
|---|---|---|---|---|---|---|
| PostOp - PreOp (no preload) | -0.5652 | -5.502 | 4.371 | 11 | -0.252 | 0.8057 |

**Shapiro-Wilk (paired differences)**: W = 0.7705, p = 0.004413 (differences deviate from normality)
