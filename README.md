# KSD on Lie Groups

This repository contains the official implementation of the experiments presented in the paper:

**Kernel Stein Discrepancy on Lie Groups: Theory and Applications**

## 1. Description

This code provides methods to compute the Minimum Kernel Stein
Discrepancy Estimator (MKSDE) and Maximum Likelihood Estimator (MLE)
of the parameter $F$ of a von-Mises Fisher (vMF) distribution defined on the
Special Orthogonal Group $SO(N)$. The experiments demonstrate the
performance of MKSDE in both parameter estimation and a
goodness-of-fit test on Lie Groups.

### Input:
- **Synthetic Data:** The code generates synthetic data on $SO(N)$ for analysis.
- **Parameter Settings:** The initial parameter values and settings are included in the scipts.

### Output: 

- **Estimates (Figure 1 in the paper):** The computed
MKSDE and MLE for the parameter $F$ and their respective distances to the ground
truth $F_0$.
- **Geodesic Distances (Table I in the paper)** The geodesic distances between the mode of estimated $F$ and the ground truth orientation of the samples in Figure 2.
- **Goodness-of-Fit Results (Table II in the paper):**
The test outputs the computed statistic $wKSD^2_n(\hat{\theta})$ and the
$(1-\beta)$-quantile.


## 2. Default Parameter Values

- **Number of Samples (n_samples):** 100 for goodness-of-fit, varies (100-200) for estimation.
- **Lie Group Dimension ($N$):** 3
- **Kernel Bandwidth ($\tau$):** 1.0
- **Initial Parameter ($F$):** 
  - $F$ for estimation: [8.5, 1.1, 4.1, 7.8, 3.9, 6, 4.3, 6.4, 4.8]
  - $F$ for goodness-of-fit: Identity matrix [1, 0, 0, 0, 1, 0, 0, 0, 1]
- **Cayley Kappa :** 2.0
- **Bootstrap Samples (n_bootstrap):** 1000
- **Significance Level ($\beta$)**: 0.1

These values can be modified within the provided scripts to suit different experimental settings.

## 3. Scripts

**Execution:** The scripts are written in R. When you execute the script, all results, including intermediate steps and computed
  values, will be displayed directly.


- `MKSDE.R`: 
	- **Functionality:** Computes the **Minimum Kernel Stein
Discrepancy Estimator (MKSDE)** and **Maximum Likelihood Estimator
(MLE)** for the parameter $F$ of a von-Mises Fisher (vMF)
distribution defined on $SO(N)$.
	- **Details:** This script estimates the
parameters $F$ using synthetic data generated from the vMF
distribution, and compares the performance of MKSDE against MLE by
computing distances and kernel values.

- `rotations/rotations*.R:`
	- **Functionality:** Computes the geodesic distances between the mode of estimated $F$ and the ground truth orientations.

	- **Details:** These scripts reproduce the results in Table II in the paper, computing the geodesic distances between the modes of MKSDE and MLE with the ground truth orientations.


- `gof.R`: 
	- **Functionality:** Performs the **MKSDE goodness-of-fit
test** to evaluate the fit of the model distribution to the data.

	- **Details:** This script generates synthetic data using the Cayley
    distribution and computes the goodness-of-fit statistics,
    including the $wKSD^2_n(\hat{\theta})$ statistic. A bootstrap procedure is used to
    obtain the $(1-\beta)$-quantile.


## Citation

If you use this codebase, or otherwise found our work valuable, please cite:

```
@article{qu2023kernel,
  title={Kernel Stein Discrepancy on Lie Groups: Theory and Applications},
  author={Qu, Xiaoda and Fan, Xiran and Vemuri, Baba C},
  journal={arXiv preprint arXiv:2305.12551},
  year={2023}
}
```

