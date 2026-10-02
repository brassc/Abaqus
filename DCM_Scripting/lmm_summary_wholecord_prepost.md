# Whole-cord linear mixed-effects model - PreOp (no preload) vs PostOp (threshold = 0.02)

Patient is a random intercept; loading condition (Flexion/Extension) is kept as its own main-effect covariate rather than averaged away - they differ hugely in magnitude, so averaging would blend two different mechanical regimes into one number. Satterthwaite-df t-tests (R `lme4`/`lmerTest`), not asymptotic z.

## Whole cord (Flexion): PreOp (no preload) vs PostOp

$$y_i = \beta_0 + \beta_{\text{statePostOp}}\,\mathbb{1}[\text{state}_i=\text{PostOp}] + u_i + \varepsilon_i$$

- $y_i$: `pct_above` (whole cord, not split by tissue) for patient $i$, Flexion only
- $\beta_0$: intercept - expected `pct_above` for PreOp (no preload) (the reference level)
- $\beta_{\text{statePostOp}}$: fixed effect of PostOp vs PreOp (no preload)
- $u_i \sim \mathcal{N}(0,\tau^2)$: random intercept per patient
- $\varepsilon_i \sim \mathcal{N}(0,\sigma^2)$: residual error

`lmer(pct_above ~ state + (1 | patient))`, Flexion data only

$$H_0:\ \beta_{\text{statePostOp}} = 0$$

No difference in % of whole cord volume above MPS $=0.02$ between PreOp (no preload) and PostOp, within Flexion.

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 17.77 | 3.34 | 22 | 5.32 | 2.436e-05 |
| statePostOp | -9.713 | 4.724 | 22 | -2.056 | 0.05182 |

## Whole cord (Extension): PreOp (no preload) vs PostOp

$$y_i = \beta_0 + \beta_{\text{statePostOp}}\,\mathbb{1}[\text{state}_i=\text{PostOp}] + u_i + \varepsilon_i$$

- $y_i$: `pct_above` (whole cord, not split by tissue) for patient $i$, Extension only
- $\beta_0$: intercept - expected `pct_above` for PreOp (no preload) (the reference level)
- $\beta_{\text{statePostOp}}$: fixed effect of PostOp vs PreOp (no preload)
- $u_i \sim \mathcal{N}(0,\tau^2)$: random intercept per patient
- $\varepsilon_i \sim \mathcal{N}(0,\sigma^2)$: residual error

`lmer(pct_above ~ state + (1 | patient))`, Extension data only

$$H_0:\ \beta_{\text{statePostOp}} = 0$$

No difference in % of whole cord volume above MPS $=0.02$ between PreOp (no preload) and PostOp, within Extension.

| Term | Estimate | Std. Error | df | t value | Pr(>\|t\|) |
|---|---|---|---|---|---|
| (Intercept) | 3.641 | 1.47 | 22 | 2.476 | 0.02144 |
| statePostOp | -0.5652 | 2.079 | 22 | -0.2718 | 0.7883 |
