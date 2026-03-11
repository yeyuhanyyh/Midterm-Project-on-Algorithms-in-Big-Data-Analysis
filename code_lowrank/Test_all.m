clear
% Generate data, choose one and uncomment it
%m = 40; n = 40; sr = 0.3; p = round(m*n*sr); r = 3; 
%m = 100; n = m; sr = 0.3; p = round(m*n*sr); r = 2;
% m = 100; n = m ; p = 5666; sr = p/m/n; r = 10;
%m = 200; n = m; p = 15665; sr = p/m/n; r = 10;
%m = 500; n = m; p = 49471; sr = p/m/n; r = 10;
m = 150; n = 300; sr = 0.49; p = round(m*n*sr); r = 10;

% Get the problem
Omega = randperm(m*n); Omega = Omega(1:p); % Omega gives the positions of the samples
xl = randn(m,r); xr = randn(n,r); A = xl*xr'; % A is the matrix to be completed
M = zeros(m, n);
M(Omega) = A(Omega); 

% Solve using approximate point gradient method and ADMM
Test_lowrank_prox;
Test_lowrank_ADMM;
