% Parameters for Augmented Lagrangian method
opts = struct('rho', 1.0, 'tol_in', 1e-10, 'tol_out', 1e-20, 'max_iter_in', 1, 'max_iter_out', 5000, 'threshold', 1e-4, 'pen', 3, 'tau', 1e3);

% Initial solution
x0 = zeros(size(A, 2), 1);

% Call Augmented Lagrangian method function
[x_alm, history] = L1ALM(x0, A, b, mu, opts);

% Compare results
fprintf('Norm of difference between CVX and ALM solutions: %e\n', norm(u - x_alm));
fprintf('Objective value for CVX solution: %e\n', cvx_optval);
fprintf('Objective value for ALM solution: %e\n', history.objval(end));

% Compute sparsity of the solution
num_nonzeros = nnz(x_alm);
total_elements = numel(x_alm);
sparsity = num_nonzeros / total_elements;

% Plot the solution
figure;
plot(x_alm);
title('Solution');

% Plot objective value history
figure;
plot(10:opts.max_iter_out, history.objval(10:end));
xlabel('Iteration');
ylabel('Objective Value');
title('Convergence of Augmented Lagrangian Method');

% Display sparsity of the solution
fprintf('Sparsity of the solution: %.2f%%\n', sparsity * 100);
