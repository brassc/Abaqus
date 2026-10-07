# PreOp (no preload) vs PostOp, GM and WM tested separately (threshold = 0.03)

GM and WM tested completely separately (does EITHER tissue's strain burden change with surgery - not whether surgery affects them differently). Each test tries a mixed model first (patient random intercept); falls back to a paired t-test if that random intercept is singular - which model actually ran is stated explicitly under each condition below. Flexion/Extension kept as separate tests, not pooled.

# Grey Matter (GM)

## GM (Flexion): PreOp (no preload) vs PostOp

**Model used: paired t-test** (LMM random-intercept variance was singular - same situation as the whole-cord comparison - so it adds nothing over a plain paired comparison).

$$H_0:\ \mu_{\Delta} = 0, \quad \Delta_i = \text{pct\_above}_{i,\text{PostOp}} - \text{pct\_above}_{i,\text{PreOp}}$$

No difference in % of GM volume above MPS $=0.03$ between PreOp (no preload) and PostOp, within Flexion (not pooled with Extension).

| Term | Estimate | CI_low | CI_high | df | t | Pr(>\|t\|) |
|---|---|---|---|---|---|---|
| PostOp - PreOp (no preload) | -2.369 | -7.213 | 2.475 | 11 | -1.076 | 0.3048 |

**Shapiro-Wilk (paired differences)**: W = 0.8408, p = 0.02831 (deviates from normality)

## GM (Extension): PreOp (no preload) vs PostOp

**Model used: linear mixed-effects model** (random intercept not singular).

$$y_i = \beta_0 + \beta_{\text{statePostOp}}\,\mathbb{1}[\text{state}_i=\text{PostOp}] + u_i + \varepsilon_i$$

$$H_0:\ \beta_{\text{statePostOp}} = 0$$

No difference in % of GM volume above MPS $=0.03$ between PreOp (no preload) and PostOp, within Extension (not pooled with Flexion).

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 0.06327 | 0.04068 | 20 | 1.555 | 0.1356 |
| statePostOp | -0.06327 | 0.05753 | 19.99 | -1.1 | 0.2845 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 8.52e-10 | 2.919e-05 |
| Residual | 0.01821 | 0.1349 |

**ICC = 0.000** (singular/near-zero)

**Shapiro-Wilk (residuals)**: W = 0.3947, p = 1.55e-08 (deviates from normality)

# White Matter (WM)

## WM (Flexion): PreOp (no preload) vs PostOp

**Model used: paired t-test** (LMM random-intercept variance was singular - same situation as the whole-cord comparison - so it adds nothing over a plain paired comparison).

$$H_0:\ \mu_{\Delta} = 0, \quad \Delta_i = \text{pct\_above}_{i,\text{PostOp}} - \text{pct\_above}_{i,\text{PreOp}}$$

No difference in % of WM volume above MPS $=0.03$ between PreOp (no preload) and PostOp, within Flexion (not pooled with Extension).

| Term | Estimate | CI_low | CI_high | df | t | Pr(>\|t\|) |
|---|---|---|---|---|---|---|
| PostOp - PreOp (no preload) | -8.589 | -17.94 | 0.7605 | 11 | -2.022 | 0.06819 |

**Shapiro-Wilk (paired differences)**: W = 0.9318, p = 0.3999

## WM (Extension): PreOp (no preload) vs PostOp

**Model used: linear mixed-effects model** (random intercept not singular).

$$y_i = \beta_0 + \beta_{\text{statePostOp}}\,\mathbb{1}[\text{state}_i=\text{PostOp}] + u_i + \varepsilon_i$$

$$H_0:\ \beta_{\text{statePostOp}} = 0$$

No difference in % of WM volume above MPS $=0.03$ between PreOp (no preload) and PostOp, within Extension (not pooled with Flexion).

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 3.129 | 0.5933 | 19.97 | 5.273 | 3.693e-05 |
| statePostOp | -2.784 | 0.8236 | 10 | -3.38 | 0.007006 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 0.1411 | 0.3757 |
| Residual | 3.731 | 1.931 |

**ICC = 0.036**

**Shapiro-Wilk (residuals)**: W = 0.9389, p = 0.188
