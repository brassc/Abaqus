### Oscillation vs no-oscillation - effect on MPS

Patient is a random intercept; thresholds are fit as separate models, not pooled. Oscillation has no Flexion/Extension split (it is its own single loading mode, not crossed with condition), unlike every other LMM in this codebase - so there is no loading-condition covariate here. Satterthwaite-df t-tests (R `lme4`/`lmerTest`), not asymptotic z.

| State | Threshold | No oscillation (mean) | Oscillation (mean) | Delta | Model | p-value | Shapiro-Wilk p |
|---|---|---|---|---|---|---|---|
| PreOp | 0.1 | 14.15 | 14.84 | +0.69 | LMM | 0.003 | 0.462 |
| PreOp | 0.15 | 7.11 | 7.57 | +0.46 | LMM | 0.043 | 0.001 |

### PreOp

#### Threshold 0.1

**Model used: linear mixed-effects model** (random intercept not singular).

$$y_i = \beta_0 + \beta_{\text{oscillation}}\,\mathbb{1}[\text{oscillation}_i=\text{Oscillation}] + u_i + \varepsilon_i$$

$$H_0:\ \beta_{\text{oscillation}} = 0$$

No difference in % of cord volume that exceeds MPS $=0.1$ at any point during the oscillation step, compared to the pre-oscillation baseline (frame 0), for PreOp.

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 14.15 | 3.89 | 11.01 | 3.637 | 0.003903 |
| oscillationOscillation | 0.694 | 0.1845 | 11 | 3.761 | 0.003151 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 181.4 | 13.47 |
| Residual | 0.2043 | 0.452 |

**ICC = 0.999**

**Shapiro-Wilk (residuals)**: W = 0.9612, p = 0.4625

#### Threshold 0.15

**Model used: linear mixed-effects model** (random intercept not singular).

$$y_i = \beta_0 + \beta_{\text{oscillation}}\,\mathbb{1}[\text{oscillation}_i=\text{Oscillation}] + u_i + \varepsilon_i$$

$$H_0:\ \beta_{\text{oscillation}} = 0$$

No difference in % of cord volume that exceeds MPS $=0.15$ at any point during the oscillation step, compared to the pre-oscillation baseline (frame 0), for PreOp.

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 7.107 | 3.063 | 11.02 | 2.32 | 0.04053 |
| oscillationOscillation | 0.4607 | 0.2019 | 11 | 2.282 | 0.0434 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 112.4 | 10.6 |
| Residual | 0.2445 | 0.4945 |

**ICC = 0.998**

**Shapiro-Wilk (residuals)**: W = 0.8388, p = 0.001365 (deviates from normality)
