### Whole cord - PreOp (no preload) vs PostOp (threshold = 0.01)

Tries a mixed model first (patient random intercept); falls back to a paired t-test if that random intercept is singular - which model actually ran is stated explicitly under each condition below (re-decided every run, not assumed from a past diagnostic). Flexion/Extension kept as separate tests, not pooled.

#### Whole cord (Flexion): PreOp (no preload) vs PostOp

**Model used: linear mixed-effects model** (random intercept not singular).

$$y_i = \beta_0 + \beta_{\text{statePostOp}}\,\mathbb{1}[\text{state}_i=\text{PostOp}] + u_i + \varepsilon_i$$

$$H_0:\ \beta_{\text{statePostOp}} = 0$$

No difference in % of whole cord volume above MPS $=0.01$ between PreOp (no preload) and PostOp, within Flexion (not pooled with Extension).

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 29.09 | 4.799 | 21.65 | 6.061 | 4.524e-06 |
| statePostOp | -15.45 | 6.339 | 11 | -2.438 | 0.03296 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 35.29 | 5.94 |
| Residual | 241.1 | 15.53 |

**ICC = 0.128**

**Shapiro-Wilk (residuals)**: W = 0.9386, p = 0.1518

#### Whole cord (Extension): PreOp (no preload) vs PostOp

**Model used: linear mixed-effects model** (random intercept not singular).

$$y_i = \beta_0 + \beta_{\text{statePostOp}}\,\mathbb{1}[\text{state}_i=\text{PostOp}] + u_i + \varepsilon_i$$

$$H_0:\ \beta_{\text{statePostOp}} = 0$$

No difference in % of whole cord volume above MPS $=0.01$ between PreOp (no preload) and PostOp, within Extension (not pooled with Flexion).

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 7.626 | 1.313 | 21.97 | 5.809 | 7.657e-06 |
| statePostOp | -7.01 | 1.823 | 11 | -3.845 | 0.002723 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 0.7342 | 0.8568 |
| Residual | 19.95 | 4.466 |

**ICC = 0.036**

**Shapiro-Wilk (residuals)**: W = 0.9223, p = 0.06567
