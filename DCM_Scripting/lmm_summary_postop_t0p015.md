### GM/WM linear mixed-effects model - PostOp (threshold = 0.015)

Patient is a random intercept; loading condition (Flexion/Extension) is kept as its own main-effect covariate rather than averaged away - they differ hugely in magnitude, so averaging would blend two different mechanical regimes into one number. Satterthwaite-df t-tests (R `lme4`/`lmerTest`), not asymptotic z.

#### Model 1 (Flexion): GM vs WM

$$y_i = \beta_0 + \beta_{\text{tissueWM}}\,\mathbb{1}[\text{tissue}_i=\text{WM}] + u_i + \varepsilon_i$$

- $y_i$: `pct_above` for patient $i$'s tissue observation (GM or WM), Flexion only
- $\beta_0$: intercept - expected `pct_above` for GM (the reference level)
- $\beta_{\text{tissueWM}}$: fixed effect of WM vs GM
- $u_i \sim \mathcal{N}(0,\tau^2)$: random intercept per patient
- $\varepsilon_i \sim \mathcal{N}(0,\sigma^2)$: residual error

`lmer(pct_above ~ tissue + (1 | patient))`, Flexion data only

$$H_0:\ \beta_{\text{tissueWM}} = 0$$

No difference in % of cord volume above MPS $=0.015$ between grey and white matter, within Flexion (not pooled with Extension).

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 6.376 | 3.379 | 12.13 | 1.887 | 0.08334 |
| tissueWM | 5.653 | 1.497 | 11 | 3.776 | 0.003069 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 123.6 | 11.12 |
| Residual | 13.45 | 3.668 |

**ICC = 0.902**

**Shapiro-Wilk (residuals)**: W = 0.926, p = 0.07933

#### Model 1 (Extension): GM vs WM

$$y_i = \beta_0 + \beta_{\text{tissueWM}}\,\mathbb{1}[\text{tissue}_i=\text{WM}] + u_i + \varepsilon_i$$

- $y_i$: `pct_above` for patient $i$'s tissue observation (GM or WM), Extension only
- $\beta_0$: intercept - expected `pct_above` for GM (the reference level)
- $\beta_{\text{tissueWM}}$: fixed effect of WM vs GM
- $u_i \sim \mathcal{N}(0,\tau^2)$: random intercept per patient
- $\varepsilon_i \sim \mathcal{N}(0,\sigma^2)$: residual error

`lmer(pct_above ~ tissue + (1 | patient))`, Extension data only

$$H_0:\ \beta_{\text{tissueWM}} = 0$$

No difference in % of cord volume above MPS $=0.015$ between grey and white matter, within Extension (not pooled with Flexion).

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 0.02292 | 0.1914 | 21.96 | 0.1198 | 0.9057 |
| tissueWM | 0.6161 | 0.265 | 11 | 2.325 | 0.04021 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 0.01824 | 0.135 |
| Residual | 0.4213 | 0.6491 |

**ICC = 0.041**

**Shapiro-Wilk (residuals)**: W = 0.7461, p = 4.373e-05 (residuals deviate from normality)
