function y_hat = simsystem(A,B,C,D,x0,u, fs,t)
%% Instructions:
% Implement a function that simulates the system here!
% Use the following function inputs and outputs.

% Function INPUT 
% A         System matrix A (matrix of size n x n)
% B         System matrix B (matrix of size n x m)
% C         System matrix C (matrix of size l x n)
% D         System matrix D (matrix of size l x m)
% x0        Initial state (vector of size n x one)
% u         system input (matrix of size N x m)

% Function OUTPUT
% yhat      predicted output (vector of size l x one)

%discrete time state space model, unknown sampling time (-1)
dt = 1/fs;
sys = ss(A,B,C,D,dt);

% %Estimation of the output
% Tfinal = (length(u)-1);
% t = 0:1/fs:Tfinal; %time

y_hat = lsim(sys,u,t,x0);

end