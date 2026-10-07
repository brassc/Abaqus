# GM/WM linear mixed-effects models (PreOp, threshold = 0.10)

Patient is a random intercept; loading condition averaged out beforehand (not modeled) to keep covariates minimal at N=12. Satterthwaite-df t-tests (R `lme4`/`lmerTest`), not asymptotic z. `tissue` and `region` are both 2-level factors, so their F-tests would just be t^2 - the t-tests below are the whole story.

## Model 1: GM vs WM

$$y_{ij} = \beta_0 + \beta_{\text{tissueWM}}\,\mathbb{1}[\text{tissue}_{ij}=\text{WM}] + u_i + \varepsilon_{ij}$$

- $y_{ij}$: `pct_above` for patient $i$'s tissue observation $j$ (GM or WM)
- $\beta_0$: intercept - expected `pct_above` for GM (the reference level)
- $\beta_{\text{tissueWM}}$: fixed effect of WM vs GM
- $u_i \sim \mathcal{N}(0,\tau^2)$: random intercept per patient (between-patient variability)
- $\varepsilon_{ij} \sim \mathcal{N}(0,\sigma^2)$: residual error

`lmer(pct_above ~ tissue + (1 | patient))`

$$H_0:\ \beta_{\text{tissueWM}} = 0$$

No difference in % of cord volume above MPS $=0.10$ between grey and white matter.

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 13.63 | 3.623 | 11.45 | 3.762 | 0.002936 |
| tissueWM | 5.318 | 1.021 | 11 | 5.207 | 0.0002911 |

## Model 2: GM vs WM, boundary vs interior

$$y_{ijk} = \beta_0 + \beta_{\text{tissueWM}}\,\mathbb{1}[\text{tissue}_{ijk}=\text{WM}] + \beta_{\text{regionInterior}}\,\mathbb{1}[\text{region}_{ijk}=\text{Interior}] + u_i + \varepsilon_{ijk}$$

- $y_{ijk}$: `pct_above` for patient $i$, tissue $j$, region $k$
- $\beta_0$: intercept - expected `pct_above` for GM x Boundary (the reference levels)
- $\beta_{\text{tissueWM}}$: fixed effect of WM vs GM, controlling for region
- $\beta_{\text{regionInterior}}$: fixed effect of Interior vs Boundary, controlling for tissue
- $u_i \sim \mathcal{N}(0,\tau^2)$: random intercept per patient
- $\varepsilon_{ijk} \sim \mathcal{N}(0,\sigma^2)$: residual error

`lmer(pct_above ~ tissue + region + (1 | patient))`

$$H_0:\ \beta_{\text{tissueWM}} = 0 \ \text{and} \ \beta_{\text{regionInterior}} = 0$$

No difference in % of cord volume above MPS $=0.10$ between grey and white matter, or between boundary and interior tissue.

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 13.63 | 3.639 | 11.49 | 3.745 | 0.003009 |
| tissueWM | 5.658 | 0.7674 | 33.01 | 7.374 | 1.809e-08 |
| regionInterior | -0.27 | 0.7674 | 33.01 | -0.3518 | 0.7272 |
