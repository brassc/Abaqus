# PreOp (no preload) vs PostOp, GM and WM tested separately (threshold = 0.01)

GM and WM tested completely separately (does EITHER tissue's strain burden change with surgery - not whether surgery affects them differently). Each test tries a mixed model first (patient random intercept); falls back to a paired t-test if that random intercept is singular - which model actually ran is stated explicitly under each condition below. Flexion/Extension kept as separate tests, not pooled.

# Grey Matter (GM)

## GM (Flexion): PreOp (no preload) vs PostOp

**Model used: paired t-test** (LMM random-intercept variance was singular - same situation as the whole-cord comparison - so it adds nothing over a plain paired comparison).

$$H_0:\ \mu_{\Delta} = 0, \quad \Delta_i = \text{pct\_above}_{i,\text{PostOp}} - \text{pct\_above}_{i,\text{PreOp}}$$

No difference in % of GM volume above MPS $=0.01$ between PreOp (no preload) and PostOp, within Flexion (not pooled with Extension).

| Term | Estimate | CI_low | CI_high | df | t | Pr(>\|t\|) |
|---|---|---|---|---|---|---|
| PostOp - PreOp (no preload) | -10.04 | -23.03 | 2.955 | 11 | -1.7 | 0.1171 |

**Shapiro-Wilk (paired differences)**: W = 0.9404, p = 0.5031

## GM (Extension): PreOp (no preload) vs PostOp

**Model used: paired t-test** (LMM random-intercept variance was singular - same situation as the whole-cord comparison - so it adds nothing over a plain paired comparison).

$$H_0:\ \mu_{\Delta} = 0, \quad \Delta_i = \text{pct\_above}_{i,\text{PostOp}} - \text{pct\_above}_{i,\text{PreOp}}$$

No difference in % of GM volume above MPS $=0.01$ between PreOp (no preload) and PostOp, within Extension (not pooled with Flexion).

| Term | Estimate | CI_low | CI_high | df | t | Pr(>\|t\|) |
|---|---|---|---|---|---|---|
| PostOp - PreOp (no preload) | 0.7082 | -4.744 | 6.16 | 11 | 0.2859 | 0.7803 |

**Shapiro-Wilk (paired differences)**: W = 0.7445, p = 0.002342 (deviates from normality)

# White Matter (WM)

## WM (Flexion): PreOp (no preload) vs PostOp

**Model used: linear mixed-effects model** (random intercept not singular).

$$y_i = \beta_0 + \beta_{\text{statePostOp}}\,\mathbb{1}[\text{state}_i=\text{PostOp}] + u_i + \varepsilon_i$$

$$H_0:\ \beta_{\text{statePostOp}} = 0$$

No difference in % of WM volume above MPS $=0.01$ between PreOp (no preload) and PostOp, within Flexion (not pooled with Extension).

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 33.6 | 5.191 | 21.51 | 6.472 | 1.833e-06 |
| statePostOp | -17.61 | 6.765 | 11 | -2.604 | 0.02454 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 48.82 | 6.987 |
| Residual | 274.6 | 16.57 |

**ICC = 0.151**

**Shapiro-Wilk (residuals)**: W = 0.9437, p = 0.1973

## WM (Extension): PreOp (no preload) vs PostOp

**Model used: paired t-test** (LMM random-intercept variance was singular - same situation as the whole-cord comparison - so it adds nothing over a plain paired comparison).

$$H_0:\ \mu_{\Delta} = 0, \quad \Delta_i = \text{pct\_above}_{i,\text{PostOp}} - \text{pct\_above}_{i,\text{PreOp}}$$

No difference in % of WM volume above MPS $=0.01$ between PreOp (no preload) and PostOp, within Extension (not pooled with Flexion).

| Term | Estimate | CI_low | CI_high | df | t | Pr(>\|t\|) |
|---|---|---|---|---|---|---|
| PostOp - PreOp (no preload) | -4.75 | -14.15 | 4.649 | 11 | -1.112 | 0.2897 |

**Shapiro-Wilk (paired differences)**: W = 0.8545, p = 0.04175 (deviates from normality)
