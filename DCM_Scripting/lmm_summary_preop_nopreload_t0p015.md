### GM/WM linear mixed-effects model - PreOp without preload (threshold = 0.015)

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
| (Intercept) | 13.77 | 4.502 | 11.84 | 3.058 | 0.01007 |
| tissueWM | 13.15 | 1.727 | 11 | 7.612 | 1.045e-05 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 225.3 | 15.01 |
| Residual | 17.9 | 4.231 |

**ICC = 0.926**

**Shapiro-Wilk (residuals)**: W = 0.9608, p = 0.4552

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
| (Intercept) | 0.9536 | 1.133 | 20.96 | 0.8415 | 0.4095 |
| tissueWM | 5.908 | 1.413 | 11 | 4.183 | 0.00153 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 3.435 | 1.853 |
| Residual | 11.97 | 3.46 |

**ICC = 0.223**

**Shapiro-Wilk (residuals)**: W = 0.9624, p = 0.488
