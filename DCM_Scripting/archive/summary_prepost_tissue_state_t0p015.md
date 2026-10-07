### PreOp (no preload) vs PostOp, GM and WM tested separately (threshold = 0.015)

GM and WM tested completely separately (does EITHER tissue's strain burden change with surgery - not whether surgery affects them differently). Each test tries a mixed model first (patient random intercept); falls back to a paired t-test if that random intercept is singular - which model actually ran is stated explicitly under each condition below. Flexion/Extension kept as separate tests, not pooled.

### Grey Matter (GM)

#### GM (Flexion): PreOp (no preload) vs PostOp

**Model used: paired t-test** (LMM random-intercept variance was singular - same situation as the whole-cord comparison - so it adds nothing over a plain paired comparison).

$$H_0:\ \mu_{\Delta} = 0, \quad \Delta_i = \text{pct\_above}_{i,\text{PostOp}} - \text{pct\_above}_{i,\text{PreOp}}$$

No difference in % of GM volume above MPS $=0.015$ between PreOp (no preload) and PostOp, within Flexion (not pooled with Extension).

| Term | Estimate | CI_low | CI_high | df | t | Pr(>\|t\|) |
|---|---|---|---|---|---|---|
| PostOp - PreOp (no preload) | -7.391 | -18.18 | 3.4 | 11 | -1.508 | 0.1598 |

**Shapiro-Wilk (paired differences)**: W = 0.9407, p = 0.5067

#### GM (Extension): PreOp (no preload) vs PostOp

**Model used: linear mixed-effects model** (random intercept not singular).

$$y_i = \beta_0 + \beta_{\text{statePostOp}}\,\mathbb{1}[\text{state}_i=\text{PostOp}] + u_i + \varepsilon_i$$

$$H_0:\ \beta_{\text{statePostOp}} = 0$$

No difference in % of GM volume above MPS $=0.015$ between PreOp (no preload) and PostOp, within Extension (not pooled with Flexion).

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 0.9536 | 0.2913 | 21.99 | 3.273 | 0.003475 |
| statePostOp | -0.9306 | 0.4083 | 11 | -2.279 | 0.04358 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 0.01812 | 0.1346 |
| Residual | 1 | 1 |

**ICC = 0.018**

**Shapiro-Wilk (residuals)**: W = 0.7331, p = 2.852e-05 (deviates from normality)

### White Matter (WM)

#### WM (Flexion): PreOp (no preload) vs PostOp

**Model used: linear mixed-effects model** (random intercept not singular).

$$y_i = \beta_0 + \beta_{\text{statePostOp}}\,\mathbb{1}[\text{state}_i=\text{PostOp}] + u_i + \varepsilon_i$$

$$H_0:\ \beta_{\text{statePostOp}} = 0$$

No difference in % of WM volume above MPS $=0.015$ between PreOp (no preload) and PostOp, within Flexion (not pooled with Extension).

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 26.91 | 4.521 | 21.87 | 5.954 | 5.562e-06 |
| statePostOp | -14.89 | 6.139 | 11 | -2.425 | 0.03373 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 19.07 | 4.367 |
| Residual | 226.2 | 15.04 |

**ICC = 0.078**

**Shapiro-Wilk (residuals)**: W = 0.9353, p = 0.1282

#### WM (Extension): PreOp (no preload) vs PostOp

**Model used: linear mixed-effects model** (random intercept not singular).

$$y_i = \beta_0 + \beta_{\text{statePostOp}}\,\mathbb{1}[\text{state}_i=\text{PostOp}] + u_i + \varepsilon_i$$

$$H_0:\ \beta_{\text{statePostOp}} = 0$$

No difference in % of WM volume above MPS $=0.015$ between PreOp (no preload) and PostOp, within Extension (not pooled with Flexion).

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 6.862 | 1.112 | 21.96 | 6.173 | 3.293e-06 |
| statePostOp | -6.223 | 1.537 | 11 | -4.049 | 0.001919 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 0.6582 | 0.8113 |
| Residual | 14.17 | 3.764 |

**ICC = 0.044**

**Shapiro-Wilk (residuals)**: W = 0.949, p = 0.2572
