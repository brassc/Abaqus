### Whole cord - PreOp (no preload) vs PostOp (threshold = 0.02)

Tries a mixed model first (patient random intercept); falls back to a paired t-test if that random intercept is singular - which model actually ran is stated explicitly under each condition below (re-decided every run, not assumed from a past diagnostic). Flexion/Extension kept as separate tests, not pooled.

#### Whole cord (Flexion): PreOp (no preload) vs PostOp

**Model used: paired t-test** (LMM random-intercept variance was singular - adds nothing over a plain paired comparison).

$$H_0:\ \mu_{\Delta} = 0, \quad \Delta_i = \text{pct\_above}_{i,\text{PostOp}} - \text{pct\_above}_{i,\text{PreOp}}$$

No difference in % of whole cord volume above MPS $=0.02$ between PreOp (no preload) and PostOp, within Flexion (not pooled with Extension).

| Term | Estimate | CI_low | CI_high | df | t | Pr(>\|t\|) |
|---|---|---|---|---|---|---|
| PostOp - PreOp (no preload) | -9.849 | -20.39 | 0.6945 | 11 | -2.056 | 0.06431 |

**Shapiro-Wilk (paired differences)**: W = 0.95, p = 0.6373

#### Whole cord (Extension): PreOp (no preload) vs PostOp

**Model used: linear mixed-effects model** (random intercept not singular).

$$y_i = \beta_0 + \beta_{\text{statePostOp}}\,\mathbb{1}[\text{state}_i=\text{PostOp}] + u_i + \varepsilon_i$$

$$H_0:\ \beta_{\text{statePostOp}} = 0$$

No difference in % of whole cord volume above MPS $=0.02$ between PreOp (no preload) and PostOp, within Extension (not pooled with Flexion).

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 3.641 | 0.601 | 21.95 | 6.058 | 4.306e-06 |
| statePostOp | -3.296 | 0.8288 | 11 | -3.976 | 0.002172 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 0.2139 | 0.4625 |
| Residual | 4.121 | 2.03 |

**ICC = 0.049**

**Shapiro-Wilk (residuals)**: W = 0.95, p = 0.2705
