### GM/WM linear mixed-effects model - PreOp with preload (threshold = 0.015)

Patient is a random intercept; loading condition (Flexion/Extension) is kept as its own main-effect covariate rather than averaged away - they differ hugely in magnitude, so averaging would blend two different mechanical regimes into one number. Satterthwaite-df t-tests (R `lme4`/`lmerTest`), not asymptotic z.

#### GM vs WM (Flexion)

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
| (Intercept) | 55.63 | 4.942 | 11.31 | 11.26 | 1.727e-07 |
| tissueWM | 7.855 | 1.165 | 11 | 6.74 | 3.2e-05 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 284.9 | 16.88 |
| Residual | 8.149 | 2.855 |

**ICC = 0.972**

**Shapiro-Wilk (residuals)**: W = 0.9866, p = 0.9812

#### GM vs WM (Extension)

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
| (Intercept) | 48.23 | 5.153 | 11.33 | 9.36 | 1.139e-06 |
| tissueWM | 6.485 | 1.25 | 11 | 5.19 | 0.0002993 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 309.2 | 17.59 |
| Residual | 9.371 | 3.061 |

**ICC = 0.971**

**Shapiro-Wilk (residuals)**: W = 0.9656, p = 0.5608
