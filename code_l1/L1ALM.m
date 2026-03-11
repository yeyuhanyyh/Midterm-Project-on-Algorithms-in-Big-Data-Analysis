function [x, out] = L1ADMM(x0, A, b, mu, opts)
% Extract options
rho = opts.rho; % Initial penalty factor
tol_in = opts.tol_in; % Inner loop stopping criterion
tol_out = opts.tol_out; % Outer loop stopping criterion
max_iter_in = opts.max_iter_in; % Maximum number of iterations for the inner loop
max_iter_out = opts.max_iter_out; % Maximum number of iterations for the outer loop
threshold = opts.threshold; % Minimum absolute value threshold (values below this are set to 0)
rho_min = opts.rho_min; % Lower bound for the penalty factor
rho_max = opts.rho_max; % Upper bound for the penalty factor

Ax = 0;
[m, n] = size(A);
x = x0; % Use the initial solution
z = zeros(m, 1);
u = zeros(m, 1);
mu_0 = mu;
mu = mu * 2^16;
alpha = 1.7;

for k = 1:max_iter_out
    % Adjust alpha after 3000 iterations
    if k > 3000
        alpha = 1;
    end
    
    % Decrease mu every 500 iterations, down to its original value
    if mod(k, 500) == 0
        mu = max(mu_0, mu / 16);
    end
    
    x_old = x; % Save x from previous iteration
    
    % x-update step using the proximal gradient method
    for i = 1:max_iter_in
        grad = A' * (Ax - z - b + u / rho);
        x_new = x - 0.1 / rho * grad;
        x_new = sign(x_new) .* max(0, abs(x_new) - 0.1 * mu / rho);
        
        % Check convergence of x-update
        if norm(x_new - x) < tol_in
            break;
        end
        
        x = x_new;
        x(abs(x) < threshold) = 0;
        Ax = A * x;
        
        % z-update step using the subgradient method
        z_new = subgradient_method(A, b, x, u, z, rho, 0.1 * mu / rho, Ax);
        
        % Save old z and check convergence of z-update
        z_old = z;
        if norm(z_new - z) < tol_in
            break;
        end
        
        z = z_new;
        % Apply relaxation with factor alpha
        z = alpha * z + (1 - alpha) * z_old;
    end
    
    % u-update step (dual variable update)
    u = u + 0.1 * (Ax - z - b) * rho;

    % Compute residuals
    r = norm(Ax - z - b);
    % s = norm(A * (x - x_old)); % Unused residual

    % Update rho if necessary
    rho = min(2 * rho, rho_max);

    % Record the history of objective value and residuals (commented out)
    % out.objval(k) = objective(A, b, mu_0, x);
    % history.r_norm(k) = r;
    % history.s_norm(k) = s;
    % out.rho(k) = rho;

    % Check convergence of outer loop
    if r < tol_out
        break;
    end
end

out.numiters = k;
end

function obj = objective(A, b, mu, x)
% Compute the objective function value
obj = mu * norm(x, 1) + norm(A * x - b, Inf);
end

function z = subgradient_method(A, b, x, u, z, rho, alpha, Ax)
% Initialization
[m, n] = size(A);

% Compute the gradient of the quadratic term
g = rho * (-Ax + b + z) - u;

% Compute the subgradient of the infinity norm
[~, I] = max(abs(z));
h = zeros(m, 1);
h(I) = sign(z(I)) / rho;

% Update z using the combined gradient and subgradient
z = -alpha / norm(g + h) * (g + h);
end

function p = proj_l1(z)
% Projection onto the L1-ball
if norm(z, 1) <= 1
    p = z;
else
    z_abs = abs(z);
    [z_sort, idx] = sort(z_abs, 'descend');
    cumulative_sum = cumsum(z_sort);
    t = (cumulative_sum - 1) ./ (1:length(z))';
    [~, i] = max(z_sort - t);
    tau = t(i);
    p = sign(z) .* max(z_abs - tau, 0);
end
end
