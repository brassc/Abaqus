### PreOp (no preload) vs PostOp, GM and WM tested separately (threshold = 0.02)

GM and WM tested completely separately (does EITHER tissue's strain burden change with surgery - not whether surgery affects them differently). Each test tries a mixed model first (patient random intercept); falls back to a paired t-test if that random intercept is singular - which model actually ran is stated explicitly under each condition below. Flexion/Extension kept as separate tests, not pooled.

### Grey Matter (GM)

#### GM (Flexion): PreOp (no preload) vs PostOp

**Model used: paired t-test** (LMM random-intercept variance was singular - same situation as the whole-cord comparison - so it adds nothing over a plain paired comparison).

$$H_0:\ \mu_{\Delta} = 0, \quad \Delta_i = \text{pct\_above}_{i,\text{PostOp}} - \text{pct\_above}_{i,\text{PreOp}}$$

No difference in % of GM volume above MPS $=0.02$ between PreOp (no preload) and PostOp, within Flexion (not pooled with Extension).

| Term | Estimate | CI_low | CI_high | df | t | Pr(>\|t\|) |
|---|---|---|---|---|---|---|
| PostOp - PreOp (no preload) | -5.309 | -13.66 | 3.038 | 11 | -1.4 | 0.1891 |

**Shapiro-Wilk (paired differences)**: W = 0.9069, p = 0.1947

#### GM (Extension): PreOp (no preload) vs PostOp

**Model used: linear mixed-effects model** (random intercept not singular).

$$y_i = \beta_0 + \beta_{\text{statePostOp}}\,\mathbb{1}[\text{state}_i=\text{PostOp}] + u_i + \varepsilon_i$$

$$H_0:\ \beta_{\text{statePostOp}} = 0$$

No difference in % of GM volume above MPS $=0.02$ between PreOp (no preload) and PostOp, within Extension (not pooled with Flexion).

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 0.4186 | 0.1638 | 21.99 | 2.555 | 0.01805 |
| statePostOp | -0.4067 | 0.2291 | 11 | -1.776 | 0.1034 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 0.007312 | 0.08551 |
| Residual | 0.3148 | 0.5611 |

**ICC = 0.023**

**Shapiro-Wilk (residuals)**: W = 0.6083, p = 7.655e-07 (deviates from normality)

### White Matter (WM)

#### WM (Flexion): PreOp (no preload) vs PostOp

**Model used: linear mixed-effects model** (random intercept not singular).

$$y_i = \beta_0 + \beta_{\text{statePostOp}}\,\mathbb{1}[\text{state}_i=\text{PostOp}] + u_i + \varepsilon_i$$

$$H_0:\ \beta_{\text{statePostOp}} = 0$$

No difference in % of WM volume above MPS $=0.02$ between PreOp (no preload) and PostOp, within Flexion (not pooled with Extension).

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 21.99 | 3.925 | 22 | 5.603 | 1.24e-05 |
| statePostOp | -12.26 | 5.52 | 11 | -2.221 | 0.04825 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 2.031 | 1.425 |
| Residual | 182.8 | 13.52 |

**ICC = 0.011**

**Shapiro-Wilk (residuals)**: W = 0.9297, p = 0.09622

#### WM (Extension): PreOp (no preload) vs PostOp

**Model used: linear mixed-effects model** (random intercept not singular).

$$y_i = \beta_0 + \beta_{\text{statePostOp}}\,\mathbb{1}[\text{state}_i=\text{PostOp}] + u_i + \varepsilon_i$$

$$H_0:\ \beta_{\text{statePostOp}} = 0$$

No difference in % of WM volume above MPS $=0.02$ between PreOp (no preload) and PostOp, within Extension (not pooled with Flexion).

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 5.091 | 0.826 | 21.95 | 6.163 | 3.37e-06 |
| statePostOp | -4.578 | 1.14 | 11 | -4.017 | 0.002025 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 0.3959 | 0.6292 |
| Residual | 7.791 | 2.791 |

**ICC = 0.048**

**Shapiro-Wilk (residuals)**: W = 0.954, p = 0.3299
