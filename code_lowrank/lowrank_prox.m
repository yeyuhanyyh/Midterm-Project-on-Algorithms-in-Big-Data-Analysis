function [X, iter, obj] = lowrank_prox(M, Omega, mu, maxIter, stepSize, tol)
% This function implements a proximal gradient method for low-rank optimization problem.
%
% Inputs:
% M: input matrix
% Omega: index set of observed entries
% mu: parameter controlling the nuclear norm of X
% maxIter: maximum number of iterations
% stepSize: step size for the gradient descent
% tol: tolerance for the stopping criterion
%
% Outputs:
% X: the solution
% iter: the number of iterations
% obj: the objective function value at each iteration

[m, n] = size(M);
mu_0=mu;
if mu==0.1
    mu=mu*2^12;
else
    mu=mu*2^16;
end

% Initialize variables
X = zeros(m,n);
% Initialize X_2 to store the last two X matrices
X_2 = zeros(m, n, 2);
obj = zeros(maxIter, 1);

% Main loop
for iter = 1:maxIter
    % Compute the gradient
    if mod(iter, 100)==0
        mu=max(mu/(2^4),mu_0);
    end
    if iter > 2
        X = X_2(:,:,2) + (iter-2)/(iter+1) * (X_2(:,:,2) - X_2(:,:,1));
    end
    gradient = 2 * P_times_X_minus_M(M, Omega, X);
    
    % Update X using gradient descent
    X = X - stepSize * gradient;
    
    % Perform singular value thresholding
    [U, S, V] = svd(X, 'econ');
    S = max(S - mu * stepSize, 0);
    X = U * S * V';
    
    % Update X_2
    if iter > 1
        X_2(:,:,1) = X_2(:,:,2);
    end
    X_2(:,:,2) = X;

    % Compute the objective value
    obj(iter) = compute_objective(M, Omega, X, mu_0);
    
    % Check convergence
    if mu==0.1
        if  iter > 300 && abs(obj(iter) - obj(iter-1))/abs(obj(iter)) < tol
            break;
        end
    else
            if iter > 400 && abs(obj(iter) - obj(iter-1))/abs(obj(iter)) < tol
        break;
            end
    end
end
end

% Helper function to compute P(X - M)
function grad = P_times_X_minus_M(M, Omega, X)
grad = zeros(size(X));
grad(Omega) = X(Omega) - M(Omega);
end

% Helper function to compute the objective value
function obj = compute_objective(M, Omega, X, mu)
obj = mu * sum(svd(X)) + norm(P_times_X_minus_M(M, Omega, X), 'fro')^2;
end
