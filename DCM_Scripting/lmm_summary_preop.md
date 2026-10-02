# GM/WM linear mixed-effects models - PreOp with preload (threshold = 0.10)

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

No difference in % of cord volume above MPS $=0.10$ between grey and white matter, within Flexion (not pooled with Extension).

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 14.21 | 3.605 | 11.8 | 3.942 | 0.002018 |
| tissueWM | 6.986 | 1.354 | 11 | 5.159 | 0.0003138 |

## Model 1 (Extension): GM vs WM

$$y_i = \beta_0 + \beta_{\text{tissueWM}}\,\mathbb{1}[\text{tissue}_i=\text{WM}] + u_i + \varepsilon_i$$

- $y_i$: `pct_above` for patient $i$'s tissue observation (GM or WM), Extension only
- $\beta_0$: intercept - expected `pct_above` for GM (the reference level)
- $\beta_{\text{tissueWM}}$: fixed effect of WM vs GM
- $u_i \sim \mathcal{N}(0,\tau^2)$: random intercept per patient
- $\varepsilon_i \sim \mathcal{N}(0,\sigma^2)$: residual error

`lmer(pct_above ~ tissue + (1 | patient))`, Extension data only

$$H_0:\ \beta_{\text{tissueWM}} = 0$$

No difference in % of cord volume above MPS $=0.10$ between grey and white matter, within Extension (not pooled with Flexion).

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 13.05 | 3.713 | 11.3 | 3.515 | 0.004649 |
| tissueWM | 3.649 | 0.8668 | 11 | 4.21 | 0.00146 |

## Model 2 (Flexion): GM vs WM, boundary vs interior

$$y_{ijk} = \beta_0 + \beta_{\text{tissueWM}}\,\mathbb{1}[\text{tissue}_{ijk}=\text{WM}] + \beta_{\text{regionInterior}}\,\mathbb{1}[\text{region}_{ijk}=\text{Interior}] + u_i + \varepsilon_{ijk}$$

- $y_{ijk}$: `pct_above` for patient $i$, tissue $j$, region $k$, Flexion only
- $\beta_0$: intercept - expected `pct_above` for GM x Boundary (the reference levels)
- $\beta_{\text{tissueWM}}$: fixed effect of WM vs GM, controlling for region
- $\beta_{\text{regionInterior}}$: fixed effect of Interior vs Boundary, controlling for tissue
- $u_i \sim \mathcal{N}(0,\tau^2)$: random intercept per patient
- $\varepsilon_{ijk} \sim \mathcal{N}(0,\sigma^2)$: residual error

`lmer(pct_above ~ tissue + region + (1 | patient))`, Flexion data only

$$H_0:\ \beta_{\text{tissueWM}} = 0 \ \text{and} \ \beta_{\text{regionInterior}} = 0$$

No difference in % of cord volume above MPS $=0.10$ between grey and white matter, or between boundary and interior tissue, within Flexion. Both factors are 2-level, so an F-test here would just be t^2 - redundant with the t-tests below.

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 14.31 | 3.622 | 11.81 | 3.952 | 0.001981 |
| tissueWM | 7.364 | 0.9685 | 33.02 | 7.604 | 9.443e-09 |
| regionInterior | -0.5217 | 0.9685 | 33.02 | -0.5387 | 0.5937 |

## Model 2 (Extension): GM vs WM, boundary vs interior

$$y_{ijk} = \beta_0 + \beta_{\text{tissueWM}}\,\mathbb{1}[\text{tissue}_{ijk}=\text{WM}] + \beta_{\text{regionInterior}}\,\mathbb{1}[\text{region}_{ijk}=\text{Interior}] + u_i + \varepsilon_{ijk}$$

- $y_{ijk}$: `pct_above` for patient $i$, tissue $j$, region $k$, Extension only
- $\beta_0$: intercept - expected `pct_above` for GM x Boundary (the reference levels)
- $\beta_{\text{tissueWM}}$: fixed effect of WM vs GM, controlling for region
- $\beta_{\text{regionInterior}}$: fixed effect of Interior vs Boundary, controlling for tissue
- $u_i \sim \mathcal{N}(0,\tau^2)$: random intercept per patient
- $\varepsilon_{ijk} \sim \mathcal{N}(0,\sigma^2)$: residual error

`lmer(pct_above ~ tissue + region + (1 | patient))`, Extension data only

$$H_0:\ \beta_{\text{tissueWM}} = 0 \ \text{and} \ \beta_{\text{regionInterior}} = 0$$

No difference in % of cord volume above MPS $=0.10$ between grey and white matter, or between boundary and interior tissue, within Extension. Both factors are 2-level, so an F-test here would just be t^2 - redundant with the t-tests below.

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 12.93 | 3.724 | 11.36 | 3.474 | 0.004971 |
| tissueWM | 3.955 | 0.6742 | 33.01 | 5.865 | 1.433e-06 |
| regionInterior | -0.01646 | 0.6742 | 33.01 | -0.02441 | 0.9807 |
