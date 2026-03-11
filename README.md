# Midterm Project on Algorithms in Big Data Analysis

This repository contains MATLAB implementations for two optimization problems from a midterm project in algorithms for big data analysis.

## Overview

The project is organized into two main parts:

1. **L1-regularized optimization / sparse recovery**
   - objective: minimize `||A*x - b||_inf + mu*||x||_1`
   - implemented with custom iterative methods and compared with solver-based baselines

2. **Low-rank matrix recovery / matrix completion**
   - recover a low-rank matrix from partially observed entries
   - implemented with proximal-gradient and ADMM-based methods

The repository also includes the final report and the original assignment handout.

## Repository Structure

```text
Midterm-Project-on-Algorithms-in-Big-Data-Analysis/
├── code_l1/
│   ├── L1ADMM.m
│   ├── L1ALM.m
│   ├── L1testdata.m
│   ├── solve_with_gurobi.m
│   ├── solve_with_mosek.m
│   ├── testL1ADMM.m
│   ├── testL1ALM.m
│   ├── testL1all.m
│   └── testL1solver.m
├── code_lowrank/
│   ├── lowrank_ADMM.m
│   ├── lowrank_prox.m
│   ├── Test_all.m
│   ├── Test_lowrank_ADMM.m
│   ├── Test_lowrank_prox.m
│   └── Test_MC.m
├── Report.pdf
├── homework-mid-req.pdf
└── README.md
```

## Part I: L1-Regularized Optimization

This part studies the optimization problem

`min_x ||A*x - b||_inf + mu*||x||_1`

which is a sparse optimization model combining an infinity-norm fitting term with an L1 regularizer.

### Implemented Methods

- `L1ADMM.m` — ADMM-based solver
- `L1ALM.m` — Augmented Lagrangian based solver
- `solve_with_mosek.m` — CVX + MOSEK baseline
- `solve_with_gurobi.m` — CVX + Gurobi baseline

### Main Scripts

- `testL1all.m` — main entry point for reproducing the full comparison
- `L1testdata.m` — synthetic data generation
- `testL1ADMM.m` — test script for the ADMM solver
- `testL1ALM.m` — test script for the ALM solver
- `testL1solver.m` — test script for the CVX baselines

### What the Script Reports

Typical outputs include:

- residual values
- L1 norms
- CPU time
- error relative to the CVX baseline

## Part II: Low-Rank Matrix Recovery

This part studies low-rank matrix recovery from partially observed entries using nuclear-norm regularization.

### Implemented Methods

- `lowrank_prox.m` — proximal-gradient method with singular value thresholding
- `lowrank_ADMM.m` — ADMM-based solver

### Main Scripts

- `Test_all.m` — main entry point for reproducing the low-rank recovery experiments
- `Test_lowrank_prox.m` — test script for the proximal-gradient method
- `Test_lowrank_ADMM.m` — test script for the ADMM method
- `Test_MC.m` — additional experiment script

### What the Script Reports

Typical outputs include:

- objective values
- iteration counts
- Frobenius norm errors
- relative reconstruction errors

## Requirements

To run the repository, you will need:

- MATLAB
- CVX (for the solver-based baseline scripts in `code_l1`)
- MOSEK or Gurobi (optional, only needed if you want to reproduce the CVX baseline comparisons)

## How to Run

### Run Part I

Open MATLAB, move to the `code_l1` folder, and run:

```matlab
testL1all
```

### Run Part II

Open MATLAB, move to the `code_lowrank` folder, and run:

```matlab
Test_all
```

## Included Documents

- `Report.pdf` — final project report
- `homework-mid-req.pdf` — original assignment description

## Notes

- The experiments are based on synthetic data generated inside the test scripts.
- Some scripts use solver-based baselines, while others implement custom first-order or splitting methods.
- Script names use the original course-project naming style.

## Author

Yuhan Ye
