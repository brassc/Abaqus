### Whole cord - PreOp (no preload) vs PostOp (threshold = 0.03)

Tries a mixed model first (patient random intercept); falls back to a paired t-test if that random intercept is singular - which model actually ran is stated explicitly under each condition below (re-decided every run, not assumed from a past diagnostic). Flexion/Extension kept as separate tests, not pooled.

#### Whole cord (Flexion): PreOp (no preload) vs PostOp

**Model used: paired t-test** (LMM random-intercept variance was singular - adds nothing over a plain paired comparison).

$$H_0:\ \mu_{\Delta} = 0, \quad \Delta_i = \text{pct\_above}_{i,\text{PostOp}} - \text{pct\_above}_{i,\text{PreOp}}$$

No difference in % of whole cord volume above MPS $=0.03$ between PreOp (no preload) and PostOp, within Flexion (not pooled with Extension).

| Term | Estimate | CI_low | CI_high | df | t | Pr(>\|t\|) |
|---|---|---|---|---|---|---|
| PostOp - PreOp (no preload) | -6.361 | -13.76 | 1.04 | 11 | -1.892 | 0.08515 |

**Shapiro-Wilk (paired differences)**: W = 0.9566, p = 0.7337

#### Whole cord (Extension): PreOp (no preload) vs PostOp

**Model used: linear mixed-effects model** (random intercept not singular).

$$y_i = \beta_0 + \beta_{\text{statePostOp}}\,\mathbb{1}[\text{state}_i=\text{PostOp}] + u_i + \varepsilon_i$$

$$H_0:\ \beta_{\text{statePostOp}} = 0$$

No difference in % of whole cord volume above MPS $=0.03$ between PreOp (no preload) and PostOp, within Extension (not pooled with Flexion).

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 2.289 | 0.3922 | 21.96 | 5.836 | 7.197e-06 |
| statePostOp | -2.055 | 0.5431 | 11 | -3.784 | 0.003027 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 0.07654 | 0.2767 |
| Residual | 1.77 | 1.33 |

**ICC = 0.041**

**Shapiro-Wilk (residuals)**: W = 0.9461, p = 0.2223
