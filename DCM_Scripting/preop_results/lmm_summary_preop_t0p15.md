### GM/WM linear mixed-effects model - PreOp with preload (threshold = 0.15)

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

No difference in % of cord volume above MPS $=0.15$ between grey and white matter, within Flexion (not pooled with Extension).

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 6.761 | 2.927 | 11.76 | 2.31 | 0.03987 |
| tissueWM | 3.917 | 1.073 | 11 | 3.651 | 0.003812 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 95.87 | 9.792 |
| Residual | 6.905 | 2.628 |

**ICC = 0.933**

**Shapiro-Wilk (residuals)**: W = 0.9801, p = 0.8979

#### GM vs WM (Extension)

$$y_i = \beta_0 + \beta_{\text{tissueWM}}\,\mathbb{1}[\text{tissue}_i=\text{WM}] + u_i + \varepsilon_i$$

- $y_i$: `pct_above` for patient $i$'s tissue observation (GM or WM), Extension only
- $\beta_0$: intercept - expected `pct_above` for GM (the reference level)
- $\beta_{\text{tissueWM}}$: fixed effect of WM vs GM
- $u_i \sim \mathcal{N}(0,\tau^2)$: random intercept per patient
- $\varepsilon_i \sim \mathcal{N}(0,\sigma^2)$: residual error

`lmer(pct_above ~ tissue + (1 | patient))`, Extension data only

$$H_0:\ \beta_{\text{tissueWM}} = 0$$

No difference in % of cord volume above MPS $=0.15$ between grey and white matter, within Extension (not pooled with Flexion).

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 6.557 | 3.009 | 11.49 | 2.179 | 0.05093 |
| tissueWM | 1.631 | 0.8916 | 11 | 1.829 | 0.09454 |

**Random effects**

| Group | Variance | Std.Dev. |
|---|---|---|
| patient | 103.9 | 10.19 |
| Residual | 4.77 | 2.184 |

**ICC = 0.956**

**Shapiro-Wilk (residuals)**: W = 0.9593, p = 0.4237
