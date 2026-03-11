% Parameters for ADMM
opts = struct('rho', 1.0, 'tol_in', 1e-10, 'tol_out', 1e-20, 'max_iter_in', 1, 'max_iter_out', 5000, 'threshold', 1e-4, 'rho_min', 0.01, 'rho_max', 1e3);

% Initial solution
x0 = zeros(size(A, 2), 1);

% Call ADMM function
[x_admm, history] = L1ADMM(x0, A, b, mu, opts);

% Compare results
fprintf('Norm of difference between true and ADMM solutions: %e\n', norm(u - x_admm));
fprintf('Objective value for CVX solution: %e\n', cvx_optval);
fprintf('Objective value for ADMM solution: %e\n', history.objval(end));

% Plot objective value history
figure;
plot(10:opts.max_iter_out, history.objval(10:end));
xlabel('Iteration');
ylabel('Objective Value');
title('Convergence of ADMM');

figure;
plot(1:opts.max_iter_out, history.rho(1:end));
xlabel('Iteration');
ylabel('rho Value');
title('Penalty');

% Compute sparsity of the solution
num_nonzeros = nnz(x_admm);
total_elements = numel(x_admm);
sparsity = num_nonzeros / total_elements;

% Display sparsity of the solution
fprintf('Sparsity of the solution: %.2f%%\n', sparsity * 100);

% Plot the solution
figure;
plot(x_admm);
title('Solution');
