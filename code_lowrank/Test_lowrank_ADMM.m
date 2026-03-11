M=zeros(m,n);
M (Omega) = A(Omega); 
norm_A = norm(A, 'fro');
X=zeros(m,n,3);
E=zeros(m,n,3);
iteration=zeros(1,3);
maxIter=100000;
StepSize=[1.5,1.5,1.5];
objective_values=cell(3);
tol=1e-6;
for index=1:size(X,3)
    para_mu =10^(-index) ;
    [X(:,:,index), iteration(index),objective_values{index}]=lowrank_ADMM(M,Omega,para_mu,1,StepSize(index),maxIter,tol);
end
% figure;
% plot(1:iteration(1), objective_values{1}(1:iteration(1)), 'LineWidth', 2);
% title('Convergence of Objective Function Values');
% xlabel('Iteration');
% ylabel('Objective Function Value');
% grid on;
for index = 1:size(X,3)
    fprintf('ADMM')
    optimal_value = min(objective_values{index}(1:iteration(index)));
    fprintf('The optimal value for mu= %d is: %f\n', 10^(-index), optimal_value);
fprintf('The number of iterations for mu= %d is: %d\n', 10^(-index), iteration(index));
error_Frobenius(index) = norm(A - X(:,:,index), 'fro');
    fprintf('The Frobenius norm of X-A for mu= %d is: %f\n', 10^(-index), error_Frobenius(index));
  relative_error_Frobenius(index) = error_Frobenius(index) / norm_A;
    fprintf('The relative Frobenius norm error for mu= %d is: %f\n', 10^(-index), relative_error_Frobenius(index));

end