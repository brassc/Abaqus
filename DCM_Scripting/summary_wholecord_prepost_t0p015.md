# Whole cord - PreOp (no preload) vs PostOp (threshold = 0.015)

Tries a mixed model first (patient random intercept); falls back to a paired t-test if that random intercept is singular - which model actually ran is stated explicitly under each condition below (re-decided every run, not assumed from a past diagnostic). Flexion/Extension kept as separate tests, not pooled.

## Whole cord (Flexion): PreOp (no preload) vs PostOp

**Model used: linear mixed-effects model** (random intercept not singular).

$$y_i = \beta_0 + \beta_{\text{statePostOp}}\,\mathbb{1}[\text{state}_i=\text{PostOp}] + u_i + \varepsilon_i$$

$$H_0:\ \beta_{\text{statePostOp}} = 0$$

No difference in % of whole cord volume above MPS $=0.015$ between PreOp (no preload) and PostOp, within Flexion (not pooled with Extension).

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 22.56 | 4.036 | 21.96 | 5.59 | 1.287e-05 |
| statePostOp | -12.38 | 5.592 | 11 | -2.214 | 0.04891 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 7.858 | 2.803 |
| Residual | 187.6 | 13.7 |

**ICC = 0.040**

**Shapiro-Wilk (residuals)**: W = 0.9319, p = 0.1074

## Whole cord (Extension): PreOp (no preload) vs PostOp

**Model used: linear mixed-effects model** (random intercept not singular).

$$y_i = \beta_0 + \beta_{\text{statePostOp}}\,\mathbb{1}[\text{state}_i=\text{PostOp}] + u_i + \varepsilon_i$$

$$H_0:\ \beta_{\text{statePostOp}} = 0$$

No difference in % of whole cord volume above MPS $=0.015$ between PreOp (no preload) and PostOp, within Extension (not pooled with Flexion).

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 4.748 | 0.8727 | 19.96 | 5.441 | 2.534e-05 |
| statePostOp | -4.323 | 1.206 | 10 | -3.584 | 0.004976 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 0.3783 | 0.6151 |
| Residual | 7.999 | 2.828 |

**ICC = 0.045**

**Shapiro-Wilk (residuals)**: W = 0.9494, p = 0.307
