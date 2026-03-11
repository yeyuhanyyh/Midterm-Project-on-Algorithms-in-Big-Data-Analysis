clear
% function Test_MC

% Minimize mu*||X||_* + sum_{ij in idx}

% Generate data
%m = 40; n = 40; sr = 0.3; p = round(m*n*sr); r = 3; 
%m = 100; n = m; sr = 0.3; p = round(m*n*sr); r = 2;
% m = 100; n = m ; p = 5666; sr = p/m/n; r = 10;
m = 200; n = m; p = 15665; sr = p/m/n; r = 10;
%m = 500; n = m; p = 49471; sr = p/m/n; r = 10;
% m = 150; n = 300; sr = 0.49; p = round(m*n*sr); r = 10;

% fr is the freedom of the set of rank-r matrices, maxr is the maximum rank one
% can recover with p samples, which is the max rank to keep fr < 1
fr = r*(m+n-r)/p; maxr = floor(((m+n)-sqrt((m+n)^2-4*p))/2);

rs = 2021; randn('state', rs); rand('state', rs);

% Get the problem
Omega = randperm(m*n); Omega = Omega(1:p); % Omega gives the positions of the samples
xl = randn(m,r); xr = randn(n,r); A = xl*xr'; % A is the matrix to be completed

% A can also be used as an individual value by indexing, selecting values column by column
M = zeros(m,n);
M(Omega) = A(Omega); % M is the sampled matrix from A

