# Midterm Project on Algorithms in Big Data Analysis

This repository contains MATLAB implementations for two optimization problems from a midterm project in algorithms for big data analysis:

1. **sparse regression with an $\ell_1$ regularizer and an $\ell_\infty$ data-fitting term**;
2. **low-rank matrix recovery with nuclear-norm regularization**.

## Repository Structure

```text
Midterm-Project-on-Algorithms-in-Big-Data-Analysis/
├── code_l1/
│   ├── L1ADMM.m
│   ├── L1ALM.m
│   ├── L1testdata.m
│   ├── ReadMe.txt
│   ├── solve_with_gurobi.m
│   ├── solve_with_mosek.m
│   ├── testL1ADMM.m
│   ├── testL1ALM.m
│   ├── testL1all.m
│   └── testL1solver.m
├── code_lowrank/
│   ├── ReadMe.txt
│   ├── Test_MC.m
│   ├── Test_all.m
│   ├── Test_lowrank_ADMM.m
│   ├── Test_lowrank_prox.m
│   ├── lowrank_ADMM.m
│   └── lowrank_prox.m
├── Report.pdf
├── homework-mid-req.pdf
└── README.md
```

## Problem I: Sparse Regression with $\ell_1$ Regularization

The first problem implemented in this repository is

$$
\min_{x \in \mathbb{R}^n} \; \|Ax-b\|_\infty + \mu \|x\|_1.
$$

Here:

- $A \in \mathbb{R}^{m \times n}$ is the data matrix,
- $b \in \mathbb{R}^m$ is the observation vector,
- $\mu > 0$ is the regularization parameter,
- $\|x\|_1 = \sum_{i=1}^n |x_i|$,
- $\|r\|_\infty = \max_i |r_i|$.

### Equivalent Split Formulation

To apply splitting-based methods such as ALM and ADMM, introduce an auxiliary variable $z$ satisfying

$$
z = Ax-b.
$$

Then the problem becomes

$$
\begin{aligned}
\min_{x,z} \quad & \|z\|_\infty + \mu \|x\|_1 \\
\text{s.t.} \quad & Ax - b = z .
\end{aligned}
$$

Equivalently,

$$
\begin{aligned}
\min_{x,z} \quad & \|z\|_\infty + \mu \|x\|_1 \\
\text{s.t.} \quad & Ax - z - b = 0 .
\end{aligned}
$$

### Augmented Lagrangian Form

A standard augmented Lagrangian for the split problem is

$$
\mathcal{L}_\rho(x,z,\lambda)
=
\mu \|x\|_1
+
\|z\|_\infty
+
\lambda^\top (Ax-z-b)
+
\frac{\rho}{2}\|Ax-z-b\|_2^2 .
$$

This is the model underlying the ALM/ADMM-type treatment in the code.

### Proximal / Shrinkage Step

The $x$-update involves the proximal operator of the $\ell_1$ norm, namely the soft-thresholding map

$$
\operatorname{shrink}(v,\tau)
=
\operatorname{sign}(v)\odot \max(|v|-\tau,0),
$$

where the max is taken componentwise.

### Synthetic Data Used in the Script

In `testL1all.m`, the numerical experiment uses randomly generated data of the form

$$
A \in \mathbb{R}^{m \times n}, \qquad
u \text{ sparse}, \qquad
b = Au,
$$

so that the recovered solution is compared against a sparse ground truth.

### Files for Problem I

- `L1ADMM.m` — ADMM-based solver
- `L1ALM.m` — Augmented Lagrangian based solver
- `solve_with_mosek.m` — solver-based baseline
- `solve_with_gurobi.m` — solver-based baseline
- `testL1all.m` — main script for reproducing the comparison

---

## Problem II: Low-Rank Matrix Recovery

The second part of the repository studies low-rank matrix recovery from partially observed entries.

Let $\Omega \subseteq \{1,\dots,m\}\times\{1,\dots,n\}$ be the set of observed indices. Define the sampling / projection operator $P_\Omega$ by

$$
(P_\Omega(X))_{ij}
=
\begin{cases}
X_{ij}, & (i,j)\in\Omega,\\
0, & (i,j)\notin\Omega.
\end{cases}
$$

The optimization problem implemented in the code is

$$
\min_{X \in \mathbb{R}^{m\times n}}
\;
\mu \|X\|_*
+
\|P_\Omega(X-M)\|_F^2,
$$

where

- $\|X\|_* = \sum_i \sigma_i(X)$ is the nuclear norm,
- $\|Y\|_F^2 = \sum_{i,j} Y_{ij}^2$ is the squared Frobenius norm,
- $M$ is the observed matrix.

### Equivalent Split Formulation for ADMM

The ADMM implementation uses the split form

$$
\begin{aligned}
\min_{X,E} \quad & \mu \|X\|_* + \|E\|_{F(\Omega)}^2 \\
\text{s.t.} \quad & X + E = M ,
\end{aligned}
$$

where

$$
\|E\|_{F(\Omega)}^2
=
\sum_{(i,j)\in\Omega} E_{ij}^2
=
\|P_\Omega(E)\|_F^2.
$$

This formulation separates the low-rank variable $X$ from the data-fitting error term.

### Proximal-Gradient Formulation

For the smooth term

$$
f(X) = \|P_\Omega(X-M)\|_F^2,
$$

the gradient is

$$
\nabla f(X) = 2P_\Omega(X-M).
$$

Therefore, a proximal-gradient step takes the form

$$
Y^k = X^k - t_k \nabla f(X^k)
    = X^k - 2t_k P_\Omega(X^k-M),
$$

followed by singular value thresholding:

$$
X^{k+1} = \operatorname{SVT}_{\mu t_k}(Y^k).
$$

If

$$
Y^k = U \operatorname{diag}(\sigma_1,\dots,\sigma_r) V^\top,
$$

then

$$
\operatorname{SVT}_\tau(Y^k)
=
U \operatorname{diag}\big((\sigma_1-\tau)_+,\dots,(\sigma_r-\tau)_+\big)V^\top,
$$

where

$$
(a)_+ = \max(a,0).
$$

### ADMM Update Form

Using the split constraint $X+E=M$, an augmented-Lagrangian style formulation can be written as

$$
\mathcal{L}_\beta(X,E,\Lambda)
=
\mu \|X\|_*
+
\|E\|_{F(\Omega)}^2
+
\langle \Lambda, X+E-M \rangle
+
\frac{1}{2\beta}\|X+E-M\|_F^2 .
$$

The code implements updates of the following form.

#### Update of $X$

The low-rank variable is updated by singular value thresholding:

$$
X^{k+1}
=
\operatorname{SVT}_{\mu \beta_k}
\left(M - E^k - \beta_k \Lambda^k\right).
$$

#### Update of $E$

Let

$$
T^k = M - X^{k+1} - \beta_k \Lambda^k.
$$

Then the error term is updated entrywise as

$$
E^{k+1}_{ij}
=
\begin{cases}
\dfrac{T^k_{ij}}{2\beta_k + 1}, & (i,j)\in\Omega,\\[8pt]
T^k_{ij}, & (i,j)\notin\Omega.
\end{cases}
$$

#### Update of the Dual Variable

The dual variable is updated by

$$
\Lambda^{k+1}
=
\Lambda^k
+
\frac{\tau}{\beta_k}
\left(X^{k+1}+E^{k+1}-M\right).
$$

### Synthetic Data Used in the Script

In `Test_all.m`, the low-rank ground-truth matrix is generated in factored form:

$$
A = X_L X_R^\top,
$$

where $X_L \in \mathbb{R}^{m\times r}$ and $X_R \in \mathbb{R}^{n\times r}$.

Then the observed matrix is formed by keeping only entries in $\Omega$:

$$
M_{ij}
=
\begin{cases}
A_{ij}, & (i,j)\in\Omega,\\
0, & (i,j)\notin\Omega.
\end{cases}
$$

### Files for Problem II

- `lowrank_prox.m` — proximal-gradient method
- `lowrank_ADMM.m` — ADMM-based method
- `Test_all.m` — main script for reproducing the experiments
- `Test_lowrank_prox.m` — test script for the proximal-gradient solver
- `Test_lowrank_ADMM.m` — test script for the ADMM solver

---

## How to Run

### Part I

```matlab
cd code_l1
testL1all
```

### Part II

```matlab
cd code_lowrank
Test_all
```

---

## Outputs

The scripts report quantities such as

- objective values,
- iteration counts,
- CPU time,
- relative errors,
- Frobenius norm reconstruction errors.

---

## Included Documents

- `Report.pdf` — final project report
- `homework-mid-req.pdf` — original assignment description

---

## Notes

This repository was created for a course project in MATLAB.  
The README is written to make the optimization models and numerical methods easier to understand for readers browsing the repository.

## Author

Yuhan Ye
