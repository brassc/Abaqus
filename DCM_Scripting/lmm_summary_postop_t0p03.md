# GM/WM linear mixed-effects model - PostOp (threshold = 0.03)

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

No difference in % of cord volume above MPS $=0.03$ between grey and white matter, within Flexion (not pooled with Extension).

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 1.61 | 1.73 | 16.73 | 0.9311 | 0.3651 |
| tissueWM | 5.097 | 1.62 | 11 | 3.146 | 0.009311 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 20.15 | 4.489 |
| Residual | 15.75 | 3.969 |

**ICC = 0.561**

**Shapiro-Wilk (residuals)**: W = 0.9179, p = 0.05249

## Model 1 (Extension): GM vs WM

$$y_i = \beta_0 + \beta_{\text{tissueWM}}\,\mathbb{1}[\text{tissue}_i=\text{WM}] + u_i + \varepsilon_i$$

- $y_i$: `pct_above` for patient $i$'s tissue observation (GM or WM), Extension only
- $\beta_0$: intercept - expected `pct_above` for GM (the reference level)
- $\beta_{\text{tissueWM}}$: fixed effect of WM vs GM
- $u_i \sim \mathcal{N}(0,\tau^2)$: random intercept per patient
- $\varepsilon_i \sim \mathcal{N}(0,\sigma^2)$: residual error

`lmer(pct_above ~ tissue + (1 | patient))`, Extension data only

$$H_0:\ \beta_{\text{tissueWM}} = 0$$

No difference in % of cord volume above MPS $=0.03$ between grey and white matter, within Extension (not pooled with Flexion).

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | -7.101e-17 | 0.1218 | 20 | -5.831e-16 | 1 |
| tissueWM | 0.3451 | 0.1722 | 19.83 | 2.004 | 0.05893 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 2.802e-09 | 5.294e-05 |
| Residual | 0.1632 | 0.4039 |

**ICC = 0.000** (singular/near-zero)

**Shapiro-Wilk (residuals)**: W = 0.6895, p = 1.445e-05 (residuals deviate from normality)
