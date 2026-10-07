# PreOp (no preload) vs PostOp, GM and WM tested separately (threshold = 0.02)

GM and WM tested completely separately (does EITHER tissue's strain burden change with surgery - not whether surgery affects them differently). Each test tries a mixed model first (patient random intercept); falls back to a paired t-test if that random intercept is singular - which model actually ran is stated explicitly under each condition below. Flexion/Extension kept as separate tests, not pooled.

# Grey Matter (GM)

## GM (Flexion): PreOp (no preload) vs PostOp

**Model used: paired t-test** (LMM random-intercept variance was singular - same situation as the whole-cord comparison - so it adds nothing over a plain paired comparison).

$$H_0:\ \mu_{\Delta} = 0, \quad \Delta_i = \text{pct\_above}_{i,\text{PostOp}} - \text{pct\_above}_{i,\text{PreOp}}$$

No difference in % of GM volume above MPS $=0.02$ between PreOp (no preload) and PostOp, within Flexion (not pooled with Extension).

| Term | Estimate | CI_low | CI_high | df | t | Pr(>\|t\|) |
|---|---|---|---|---|---|---|
| PostOp - PreOp (no preload) | -5.195 | -13.57 | 3.179 | 11 | -1.365 | 0.1994 |

**Shapiro-Wilk (paired differences)**: W = 0.9013, p = 0.1649

## GM (Extension): PreOp (no preload) vs PostOp

**Model used: linear mixed-effects model** (random intercept not singular).

$$y_i = \beta_0 + \beta_{\text{statePostOp}}\,\mathbb{1}[\text{state}_i=\text{PostOp}] + u_i + \varepsilon_i$$

$$H_0:\ \beta_{\text{statePostOp}} = 0$$

No difference in % of GM volume above MPS $=0.02$ between PreOp (no preload) and PostOp, within Extension (not pooled with Flexion).

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 0.4186 | 0.9932 | 21.98 | 0.4215 | 0.6775 |
| statePostOp | 1.614 | 1.386 | 11 | 1.165 | 0.2688 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 0.3143 | 0.5606 |
| Residual | 11.52 | 3.395 |

**ICC = 0.027**

**Shapiro-Wilk (residuals)**: W = 0.6026, p = 6.603e-07 (deviates from normality)

# White Matter (WM)

## WM (Flexion): PreOp (no preload) vs PostOp

**Model used: paired t-test** (LMM random-intercept variance was singular - same situation as the whole-cord comparison - so it adds nothing over a plain paired comparison).

$$H_0:\ \mu_{\Delta} = 0, \quad \Delta_i = \text{pct\_above}_{i,\text{PostOp}} - \text{pct\_above}_{i,\text{PreOp}}$$

No difference in % of WM volume above MPS $=0.02$ between PreOp (no preload) and PostOp, within Flexion (not pooled with Extension).

| Term | Estimate | CI_low | CI_high | df | t | Pr(>\|t\|) |
|---|---|---|---|---|---|---|
| PostOp - PreOp (no preload) | -12.11 | -24.28 | 0.0575 | 11 | -2.191 | 0.05092 |

**Shapiro-Wilk (paired differences)**: W = 0.9246, p = 0.3265

## WM (Extension): PreOp (no preload) vs PostOp

**Model used: paired t-test** (LMM random-intercept variance was singular - same situation as the whole-cord comparison - so it adds nothing over a plain paired comparison).

$$H_0:\ \mu_{\Delta} = 0, \quad \Delta_i = \text{pct\_above}_{i,\text{PostOp}} - \text{pct\_above}_{i,\text{PreOp}}$$

No difference in % of WM volume above MPS $=0.02$ between PreOp (no preload) and PostOp, within Extension (not pooled with Flexion).

| Term | Estimate | CI_low | CI_high | df | t | Pr(>\|t\|) |
|---|---|---|---|---|---|---|
| PostOp - PreOp (no preload) | -1.514 | -7.425 | 4.398 | 11 | -0.5636 | 0.5843 |

**Shapiro-Wilk (paired differences)**: W = 0.807, p = 0.01126 (deviates from normality)
