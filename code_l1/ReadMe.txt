solve_with_mosek is a function that uses cvx to call mosek to solve the problem directly;
solve_with_gurobi is a function that uses cvx to call gurobi to solve the problem directly;

%%%%%%%%%%%%%

L1ALM is a function that solves the original problem using the Augmented Lagrangian Method. The interface is as follows:

[x, out] = L1ALM(x0, A, b, mu, opts)

The elements of the opts structure are as follows: % opt.rho Initial penalty factor
% opt.tol_in Stopping criterion for the inner loop
% opt.tol_out Stopping criterion for the outer loop
% opt.max_iter_in Maximum number of iterations for the inner loop
% opt.max_iter_out Maximum number of iterations for the outer loop
% opt.threshold Minimum absolute value threshold (values below this are set to 0)
% opt.pen Penalty factor update coefficient
% opt.rho_max Upper bound for the penalty factor

The recommended settings are:
opts1 = struct('rho', 1.0, 'tol_in', 1e-20, 'tol_out', 1e-7, 'max_iter_in', 1, 'max_iter_out', 4000, 'threshold', 1e-4, 'pen', 2, 'rho_max', 1e3);
Here, tol_out=1e-7 ensures precision while keeping the algorithm relatively fast.

% out.numiters records the number of iterations performed by the algorithm.

%%%%%%%%%%%%%%

L1ADMM is a function that solves the original problem using the Alternating Direction Method of Multipliers (ADMM). The interface is as follows:

[x, out] = L1ALM(x0, A, b, mu, opts)

The elements of the opts structure are as follows:
% opts.rho Initial penalty factor
% opts.tol_in Stopping criterion for the inner loop
% opts.tol_out Stopping criterion for the outer loop
% opts.max_iter_in Maximum number of iterations for the inner loop
% opts.max_iter_out Maximum number of iterations for the outer loop
% opts.threshold Minimum absolute value threshold (values below this are set to 0)
% opts.rho_min Lower bound for the penalty factor
% opts.rho_max Upper bound for the penalty factor

The recommended settings are:
opts2 = struct('rho', 1.0, 'tol_in', 1e-20, 'tol_out', 1e-7, 'max_iter_in', 1, 'max_iter_out', 4000, 'threshold', 1e-4, 'rho_min', 0.01, 'rho_max', 1e3);
Here, tol_out=1e-7 ensures precision while keeping the algorithm relatively fast.

% out.numiters records the number of iterations performed by the algorithm.

%%%%%%%%%%%%%%

Run testL1all to obtain all results (time, error, etc.)