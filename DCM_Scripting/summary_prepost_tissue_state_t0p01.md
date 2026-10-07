# PreOp (no preload) vs PostOp, GM and WM tested separately (threshold = 0.01)

GM and WM tested completely separately (does EITHER tissue's strain burden change with surgery - not whether surgery affects them differently). Each test tries a mixed model first (patient random intercept); falls back to a paired t-test if that random intercept is singular - which model actually ran is stated explicitly under each condition below. Flexion/Extension kept as separate tests, not pooled.

# Grey Matter (GM)

## GM (Flexion): PreOp (no preload) vs PostOp

**Model used: linear mixed-effects model** (random intercept not singular).

$$y_i = \beta_0 + \beta_{\text{statePostOp}}\,\mathbb{1}[\text{state}_i=\text{PostOp}] + u_i + \varepsilon_i$$

$$H_0:\ \beta_{\text{statePostOp}} = 0$$

No difference in % of GM volume above MPS $=0.01$ between PreOp (no preload) and PostOp, within Flexion (not pooled with Extension).

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 19.86 | 4.191 | 22 | 4.739 | 9.944e-05 |
| statePostOp | -10.31 | 5.885 | 11 | -1.752 | 0.1076 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 2.92 | 1.709 |
| Residual | 207.8 | 14.42 |

**ICC = 0.014**

**Shapiro-Wilk (residuals)**: W = 0.921, p = 0.06142

## GM (Extension): PreOp (no preload) vs PostOp

**Model used: paired t-test** (LMM random-intercept variance was singular - same situation as the whole-cord comparison - so it adds nothing over a plain paired comparison).

$$H_0:\ \mu_{\Delta} = 0, \quad \Delta_i = \text{pct\_above}_{i,\text{PostOp}} - \text{pct\_above}_{i,\text{PreOp}}$$

No difference in % of GM volume above MPS $=0.01$ between PreOp (no preload) and PostOp, within Extension (not pooled with Flexion).

| Term | Estimate | CI_low | CI_high | df | t | Pr(>\|t\|) |
|---|---|---|---|---|---|---|
| PostOp - PreOp (no preload) | -2.315 | -4.185 | -0.4454 | 10 | -2.759 | 0.02016 |

**Shapiro-Wilk (paired differences)**: W = 0.8204, p = 0.01748 (deviates from normality)

# White Matter (WM)

## WM (Flexion): PreOp (no preload) vs PostOp

**Model used: linear mixed-effects model** (random intercept not singular).

$$y_i = \beta_0 + \beta_{\text{statePostOp}}\,\mathbb{1}[\text{state}_i=\text{PostOp}] + u_i + \varepsilon_i$$

$$H_0:\ \beta_{\text{statePostOp}} = 0$$

No difference in % of WM volume above MPS $=0.01$ between PreOp (no preload) and PostOp, within Flexion (not pooled with Extension).

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 33.6 | 5.231 | 21.4 | 6.422 | 2.099e-06 |
| statePostOp | -17.94 | 6.75 | 11 | -2.658 | 0.02226 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 55 | 7.416 |
| Residual | 273.4 | 16.53 |

**ICC = 0.167**

**Shapiro-Wilk (residuals)**: W = 0.9404, p = 0.1663

## WM (Extension): PreOp (no preload) vs PostOp

**Model used: linear mixed-effects model** (random intercept not singular).

$$y_i = \beta_0 + \beta_{\text{statePostOp}}\,\mathbb{1}[\text{state}_i=\text{PostOp}] + u_i + \varepsilon_i$$

$$H_0:\ \beta_{\text{statePostOp}} = 0$$

No difference in % of WM volume above MPS $=0.01$ between PreOp (no preload) and PostOp, within Extension (not pooled with Flexion).

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 9.399 | 1.808 | 19.98 | 5.199 | 4.372e-05 |
| statePostOp | -8.551 | 2.514 | 10 | -3.402 | 0.006749 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 1.198 | 1.094 |
| Residual | 34.75 | 5.895 |

**ICC = 0.033**

**Shapiro-Wilk (residuals)**: W = 0.9187, p = 0.07152
