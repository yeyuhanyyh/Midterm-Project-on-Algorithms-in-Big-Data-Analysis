cvx_solver mosek
cvx_begin
    variable x(n)
    minimize( norm(A*x - b, inf) + mu*norm(x,1) )
cvx_end
cvx_solver gurobi
cvx_begin
    variable x(n)
    minimize( norm(A*x - b, inf) + mu*norm(x,1) )
cvx_end
fprintf('Objective value for CVX solution: %e\n', cvx_optval);