# GM/WM linear mixed-effects model - PostOp (threshold = 0.01)

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
| (Intercept) | 9.82 | 4.199 | 11.74 | 2.339 | 0.03791 |
| tissueWM | 6.164 | 1.514 | 11 | 4.072 | 0.001844 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 197.8 | 14.06 |
| Residual | 13.75 | 3.708 |

**ICC = 0.935**

**Shapiro-Wilk (residuals)**: W = 0.8755, p = 0.006744 (residuals deviate from normality)

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
| (Intercept) | 3.355 | 2.666 | 11.55 | 1.258 | 0.2331 |
| tissueWM | 1.733 | 0.8295 | 11 | 2.089 | 0.06074 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 81.17 | 9.01 |
| Residual | 4.128 | 2.032 |

**ICC = 0.952**

**Shapiro-Wilk (residuals)**: W = 0.8106, p = 0.0004405 (residuals deviate from normality)
