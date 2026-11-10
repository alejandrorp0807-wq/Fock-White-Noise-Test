# Fock-Serial-Independence

### Detecting nonlinear and higher-order temporal dependence using coherent-state representations

**Author:** Alejandro Rodríguez Peláez  
**Status:** Exploratory research project · November 2026

## Overview

**Fock-Serial-Independence** investigates a hierarchy of statistical tests for serial independence based on finite-dimensional coherent-state representations in Fock space. Instead of examining only ordinary autocorrelations, the method measures interactions between transformed observations at two, three, and four time points.

The goal is to explore whether **modal representations** offer interpretable diagnostics for nonlinear and higher-order dependence, particularly when pairwise correlations vanish.

> **Research caveat:** This is an exploratory implementation, **not** a validated replacement for established tests such as BDS, HSIC, or Ljung–Box. The available simulations do not establish general superiority or a novel universal independence test.

## Mathematical idea

For a real-valued observation, consider the coherent state

$$
|\alpha(x)\rangle = e^{-\alpha(x)^2/2}\sum_{n=0}^{\infty}\frac{\alpha(x)^n}{\sqrt{n!}}|n\rangle.
$$

A finite expansion gives features

$$
f_n(x)=e^{-\alpha(x)^2/2}\frac{\alpha(x)^n}{\sqrt{n!}},\qquad n=0,\ldots,K.
$$

The code applies a rank-to-normal transformation, centers the modal features, and optionally scales each coordinate. It then constructs lagged interaction tensors. For example, a three-time interaction has components

$$
\widehat{T}_{ijk}=\frac{1}{N-2}\sum_{t=3}^{N}\widetilde f_i(X_t)\widetilde f_j(X_{t-1})\widetilde f_k(X_{t-2}).
$$

Tests use either a squared tensor norm or, in an **exploratory extension**, the largest absolute modal interaction. P-values are calibrated by randomly permuting the observed sequence.

In the infinite-mode limit, the inner product of coherent states reproduces a Gaussian kernel. Consequently, the method is closely connected to **kernel-based dependence testing**, rather than being an entirely independent mathematical framework.

## What is being tested?

- **Null hypothesis:** observations are independent and identically distributed (i.i.d.), or exchangeable under the permutation framework.
- **Alternative:** temporal dependence detectable through the selected modal interactions.
- **Important:** rejecting i.i.d. is **not equivalent** to rejecting *weak white noise*. An ARCH process, for example, may be uncorrelated but serially dependent.

The default primary analysis preselects **K = 4**, a bandwidth of **1.4**, pairwise lags 1–5, and a consecutive three-time interaction. The primary combined p-value uses a two-test Bonferroni adjustment. Exploring additional mode orders, bandwidths, or temporal interaction orders requires appropriate multiplicity control.

## Main experimental findings

The Monte Carlo study compares the modal tests with **Ljung–Box**, **ARCH-LM**, and **BDS** under i.i.d. Gaussian/heavy-tailed/non-Gaussian noise, AR(1), ARCH, bilinear dependence, stochastic volatility, and specially constructed three- and four-time dependencies.

With **N = 1,500** observations and **35 replications per scenario** in the main study:

| Scenario | Ljung–Box | BDS (up to order 4) | Fock, 2-time (K=4) | Fock, 3-time (K=4) |
|---|---:|---:|---:|---:|
| Gaussian i.i.d. (false-positive rate) | 14% | 3% | 6% | 6% |
| AR(1) | 100% | 100% | 100% | 89% |
| ARCH | 49% | 100% | 100% | 100% |
| Noisy three-time dependence | 3% | 100% | 3% | 100% |
| Noisy four-time dependence | 3% | 100% | 0% | 6% |

**Interpretation:** Higher-order modal interactions detect dependencies that conventional autocorrelation misses. However, **BDS performed at least as well in these benchmark scenarios**. The post-hoc maximum-interaction variant performed much better than the tensor norm for certain sparse higher-order alternatives, but this result is exploratory and also reproducible with simpler targeted sign-product tests.

All percentages are Monte Carlo estimates, not population power or formal guarantees; 35 replications entail substantial uncertainty.

## Installation

Requires **Python 3.10+** (recommended) and the dependencies in `requirements.txt`.

```bash
python -m venv .venv
# Linux/macOS:
source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

## Quick start: test your own series

```python
import numpy as np
from fock_white_noise_study import permutation_pvalues

rng = np.random.default_rng(123)
x = rng.normal(size=1_200)

results = permutation_pvalues(
    x,
    Ks=(4,),
    bandwidth=1.4,
    perms=999,
    seed=2026,
)

print("Primary adjusted p-value:", results["primary_Fock_K4"])
print("Two-time p-value:", results["fock_q2_K4"])
print("Three-time p-value:", results["fock_q3_K4"])
```

Choose the settings **before** inspecting a dataset. For a family of six Bonferroni-adjusted tests, 99 permutations are insufficient to reach adjusted p < 0.05 because the smallest possible p-value is 0.01.

## Reproduce the Monte Carlo study

```bash
# Main benchmark; results are written into dedicated folders.
python fock_white_noise_study.py --suite main --out run_main600 --sizes 600 --reps 55 --perms 99 --with-bds
python fock_white_noise_study.py --suite main --out run_main1500 --sizes 1500 --reps 35 --perms 99 --with-bds

# Four-time interaction study.
python fock_white_noise_study.py --suite quad --out run_quad --sizes 600 --reps 30 --perms 99 --with-bds

# Exploratory sparse maximum-interaction test.
python sparse_modal_followup.py --n 600 --reps 35 --perms 199 --out run_sparse

# Validation checks.
python unit_tests.py
```

For additional tests and plotting scripts, consult the accompanying manuscript and the study's reproduction notes. Computational cost increases substantially for high-order tensors and many permutations.

## Repository files

| File or directory | Description |
|---|---|
| `fock_white_noise_study.py` | Main coherent-state tests, Monte Carlo scenarios, permutation calibration, classical benchmarks |
| `sparse_modal_followup.py` | Exploratory maximum-interaction statistic |
| `higher_order_controls.py` | Targeted sign-product benchmarks |
| `unit_tests.py` | Basic numerical and construction checks |
| `empirical_demo.py` | Illustrative empirical diagnostic; not independent validation |
| `requirements.txt` | Python dependencies |
| `results/` | Simulation-level CSVs and aggregated power estimates |
| `figs/` | Figures used in the manuscript |
| `manuscrito.tex` | LaTeX manuscript source |
| `Manuscrito_Test_Independencia_Fock_2026.pdf` | Research manuscript (Spanish); an English edition is available separately |

## Scientific limitations

1. **Originality:** coherent-state inner products give Gaussian kernels, and kernel-based independence tests already exist.
2. **Power:** superiority over BDS, dHSIC/HSIC, or other modern competitors has **not** been established.
3. **Post-hoc exploration:** the sparse maximum-interaction extension was developed after inspecting initial simulation results and needs independent confirmation.
4. **Calibration:** permutation validity relies on the null's exchangeability assumptions and preselected test configurations; model-fitted residuals and nonstationary data require extra care.
5. **Computational scaling:** high temporal interaction orders grow rapidly with the number of features; cost and runtime must be compared fairly.
6. **Multiple testing:** tuning mode count, lag set, bandwidth, and interaction order on the same sample can invalidate nominal p-values without appropriate corrections.

## Manuscript and citation

**Alejandro Rodríguez Peláez (2026).** *Coherent-State Fock Representations for Detecting Nonlinear Temporal Dependencies: Theory and Monte Carlo Evaluation.* Exploratory research manuscript.

If you reuse this work, please cite the manuscript and link to the repository once it is public. Do not interpret this repository as peer-reviewed evidence of a superior statistical test.

## License

**License not yet specified.** Add a `LICENSE` file before inviting reuse or contributions.
