clear;
rng(1);
% Test data
n = 1024;
m = 512;
A = randn(m, n);
u = sprandn(n, 1, 0.1);
b = A * u;
mu = 1e-2;

% Initial solutions
x0 = rand(n, 1);

% Using mosek as cvx_solver
[x1, ~] = solve_with_mosek(x0, A, b, mu);

% Using gurobi as cvx_solver
[x2, ~] = solve_with_gurobi(x0, A, b, mu);

% Parameters for ALM
opts1 = struct('rho', 1.0, 'tol_in', 1e-20, 'tol_out', 1e-7, 'max_iter_in', 1, 'max_iter_out', 4000, 'threshold', 1e-4, 'pen', 2,'rho_min', 0.01, 'rho_max', 1e3);
% rho % Initial penalty factor
% tol_in % Inner loop stopping condition
% tol_out % Outer loop stopping condition
% max_iter_in % Maximum number of iterations for the inner loop
% max_iter_out % Maximum number of iterations for the outer loop
% threshold % Minimum absolute value threshold (values below this will be set to 0)
% pen % Penalty factor update coefficient
% rho_max % Upper bound for penalty factor

% Parameters for ADMM
opts2 = struct('rho', 1.0, 'tol_in', 1e-20, 'tol_out', 1e-7, 'max_iter_in', 1, 'max_iter_out', 4000, 'threshold', 1e-4, 'rho_min', 0.01, 'rho_max', 1e3);
% rho Initial penalty factor
% tol_in Inner loop stopping condition
% tol_out Outer loop stopping condition (verifying primal feasibility)
% max_iter_in Maximum number of iterations for the inner loop
% max_iter_out Maximum number of iterations for the outer loop
% threshold Minimum absolute value threshold (values below this will be set to 0)
% rho_min Lower bound for penalty factor
% rho_max Upper bound for penalty factor

% ADMM
tic;
[x3, out3] = L1ADMM(x0, A, b, mu, opts2);
t3 = toc;

% ALM
tic;
[x4, out4] = L1ALM(x0, A, b, mu, opts1);
t4 = toc;

% Error functions
errfun = @(x1, x2) norm(x1 - x2) / (1 + norm(x1));
resfun = @(x) norm(A * x - b, 1);
nrm1fun = @(x) norm(x, 1);

% Print comparison results
fprintf('cvx-mosek:   res: %3.2e, nrm1: %3.2e\n', resfun(x1), nrm1fun(x1));
fprintf('cvx-gurobi:  res: %3.2e, nrm1: %3.2e\n', resfun(x2), nrm1fun(x2));
fprintf('ALM:         res: %3.2e, nrm1: %3.2e, cpu: %5.2f, err-to-cvx-mosek: %3.2e\n', ...
        resfun(x4), nrm1fun(x4), t4, errfun(x1, x4));
fprintf('ADMM:        res: %3.2e, nrm1: %3.2e, cpu: %5.2f, err-to-cvx-mosek: %3.2e\n', ...
        resfun(x3), nrm1fun(x3), t3, errfun(x1, x3));


