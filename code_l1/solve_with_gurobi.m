% Function to use Gurobi as the cvx_solver
function [x, out] = solve_with_gurobi(x0, A, b, mu, opts)
    n = length(x0);
    cvx_solver gurobi
    cvx_begin quiet
        variable x(n)
        minimize(norm(A*x - b, inf) + mu*norm(x, 1))
    cvx_end
    out = cvx_optval;
end
