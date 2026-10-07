# GM/WM linear mixed-effects model - PreOp without preload (threshold = 0.01)

Patient is a random intercept; loading condition (Flexion/Extension) is kept as its own main-effect covariate rather than averaged away - they differ hugely in magnitude, so averaging would blend two different mechanical regimes into one number. Satterthwaite-df t-tests (R `lme4`/`lmerTest`), not asymptotic z.

## Model 1 (Flexion): GM vs WM

$$y_i = \beta_0 + \beta_{\text{tissueWM}}\,\mathbb{1}[\text{tissue}_i=\text{WM}] + u_i + \varepsilon_i$$

- $y_i$: `pct_above` for patient $i$'s tissue observation (GM or WM), Flexion only
- $\beta_0$: intercept - expected `pct_above` for GM (the reference level)
- $\beta_{\text{tissueWM}}$: fixed effect of WM vs GM
- $u_i \sim \mathcal{N}(0,\tau^2)$: random intercept per patient
- $\varepsilon_i \sim \mathcal{N}(0,\sigma^2)$: residual error

`lmer(pct_above ~ tissue + (1 | patient))`, Flexion data only

$$H_0:\ \beta_{\text{tissueWM}} = 0$$

No difference in % of cord volume above MPS $=0.01$ between grey and white matter, within Flexion (not pooled with Extension).

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 19.86 | 5.166 | 11.55 | 3.844 | 0.002498 |
| tissueWM | 13.74 | 1.617 | 11 | 8.499 | 3.657e-06 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 304.5 | 17.45 |
| Residual | 15.68 | 3.96 |

**ICC = 0.951**

**Shapiro-Wilk (residuals)**: W = 0.9487, p = 0.2541

## Model 1 (Extension): GM vs WM

$$y_i = \beta_0 + \beta_{\text{tissueWM}}\,\mathbb{1}[\text{tissue}_i=\text{WM}] + u_i + \varepsilon_i$$

- $y_i$: `pct_above` for patient $i$'s tissue observation (GM or WM), Extension only
- $\beta_0$: intercept - expected `pct_above` for GM (the reference level)
- $\beta_{\text{tissueWM}}$: fixed effect of WM vs GM
- $u_i \sim \mathcal{N}(0,\tau^2)$: random intercept per patient
- $\varepsilon_i \sim \mathcal{N}(0,\sigma^2)$: residual error

`lmer(pct_above ~ tissue + (1 | patient))`, Extension data only

$$H_0:\ \beta_{\text{tissueWM}} = 0$$

No difference in % of cord volume above MPS $=0.01$ between grey and white matter, within Extension (not pooled with Flexion).

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 2.647 | 1.755 | 17.77 | 1.508 | 0.1491 |
| tissueWM | 7.191 | 1.776 | 11 | 4.048 | 0.001923 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 18.02 | 4.245 |
| Residual | 18.94 | 4.351 |

**ICC = 0.488**

**Shapiro-Wilk (residuals)**: W = 0.9444, p = 0.204
