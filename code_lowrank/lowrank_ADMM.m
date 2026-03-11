function [X, iter, obj] = lowrank_ADMM(M, Omega, mu, para_beta, tau, maxIter, tol)
% This function implements the ADMM algorithm for solving the optimization problem:
% min_{X,E} mu*||X||_* + ||E||_{F(Omega)}^2, subject to X + E = M
%
% Inputs:
% M: input matrix
% Omega: index set of observed entries
% mu: parameter controlling the nuclear norm of X
% para_beta: ADMM penalty parameter
% StepSize: step size for the update of the dual variable
% maxIter: maximum number of iterations
% tol: tolerance for the stopping criterion
%
% Outputs:
% X: the solution
% iter: the number of iterations
% obj: the objective function value at each iteration

% Initialize variables
[m, n] = size(M);
X = randn(m, n);
E = randn(m, n);
Lambda = randn(m, n);
mu_0=mu;
mu=mu*2^16;
% ADMM parameters
rho = 1;  % Parameter for updating the penalty parameter

% Initialize iteration counter
iter = 0;

% Main loop
while iter < maxIter
    iter = iter + 1;
    if mod(iter, 100)==0
        mu=max(mu/(2^4),mu_0);
    end
    % Update X
    X_old = X;
    X = update_x(E, Lambda, M, mu, para_beta);
    
    % Update E
    E_new = solveE(M, X, Lambda, Omega, para_beta);
    E = E_new;
    
    % Update the dual variable
    Lambda = Lambda + tau/para_beta * (X + E - M);
    
    % Compute the objective value
    obj(iter) = compute_objective(M, Omega, X, mu_0);
    
    % Check primal feasibility
    if norm(X + E - M, 'fro') / norm(M, 'fro') < 1e-100
        break;
    end
    
    % Check convergence
    if iter > 400 && abs(obj(iter) - obj(iter-1))/abs(obj(iter)) < tol
        break;
    end
    
    % Update the penalty parameter
    para_beta = rho * para_beta;
end
end

% The following are helper functions used in the main function

function X_new = update_x(E, Lambda, M, mu, para_beta)
% This function updates X using singular value thresholding
s_j = M - E - (para_beta) * Lambda;
[U, Sigma, V] = svd(s_j, 'econ');
diagSigma = diag(Sigma);
diagSigma = shrink(diagSigma, mu*para_beta);
X_new = U * diag(diagSigma) * V';
end

function y = shrink(x, tau)
% This function implements the shrinkage operator
y = sign(x) .* max(abs(x) - tau, 0);
end

function E = solveE(M, X, Lambda, Omega, para_beta)
% This function solves the subproblem for E
[m, n] = size(M);
T = M - X - (para_beta) * Lambda;
E = T;
E(Omega) = T(Omega)./(2*para_beta+1);
end

function obj = compute_objective(M, Omega, X, mu)
% This function computes the objective function value
obj = mu * sum(svd(X)) + norm(P_times_X_minus_M(M, Omega, X), 'fro')^2;
end

function grad = P_times_X_minus_M(M, Omega, X)
% This function computes P(X - M)
grad = zeros(size(X));
grad(Omega) = X(Omega) - M(Omega);
end
