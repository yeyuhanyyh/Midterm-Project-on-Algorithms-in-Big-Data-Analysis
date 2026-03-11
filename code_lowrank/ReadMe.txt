Run Test_all to get all results (optimal values, number of iterations, etc.)
lowrank_prox is a 3.1 function that solves the low-rank matrix recovery problem using the approximate point gradient method;
The function interface is as follows:
function [X, iter, obj] = lowrank_prox(M, Omega, mu, maxIter, stepSize, tol)
% Inputs.
% M: input matrix
% Omega: index set of observed entries
% mu: parameter controlling the nuclear norm of X
% maxIter: maximum number of iterations
% stepSize: step size for the gradient descent
% tol: tolerance for the stopping criterion
See Test_lowrank_prox for details on the optimal parameters, or just run Test_all (default optimal parameters)
% Outputs: tolerance for the stopping criterion
% Outputs.
% X: the solution
% iter: the number of iterations
% obj: the objective function value at each iteration
lowrank_ADMM is a function of 3.2 that solves the low-rank matrix recovery problem using ADMM;
The function interface is as follows:
function [X, iter, obj] = lowrank_ADMM(M, Omega, mu, para_beta, tau, maxIter, tol)
% Inputs.
% M: input matrix
% Omega: index set of observed entries
% mu: parameter controlling the nuclear norm of X
% para_beta: ADMM penalty parameter
% StepSize: step size for the update of the dual variable
% maxIter: maximum number of iterations
% tol: tolerance for the stopping criterion
See Test_lowrank_ADMM for details of the optimization parameters, or just run Test_all (default optimization parameters)
% Outputs: maximum number of iterations % tol: tolerance for the stopping criterion
% Outputs.
% X: the solution
% iter: the number of iterations
% obj: the objective function value at each iteration

Translated with DeepL.com (free version)