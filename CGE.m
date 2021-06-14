%% Graph envelope estimation using convex optimization
% Arun Venkitaraman 23-03-2017
% This code computes the convex graph envelope of a signal
% Inputs: x: signal, Anorm: normalized adj matix, beta: regularization
% paramter
% Output: x_env, CGE of x, and corresponding graph phase modulation x_pm


function [x_env, x_pm]=CGE(x,Anorm,bet)
%x_env:  Convex graph envelope


N=length(x);
cvx_begin quiet
variable x_env(N,1) nonnegative;

minimize( quad_form(x_env-Anorm*x_env,eye(N))+bet*MS_g(x_env,Anorm,2));
subject to

 x_env >= abs(x); 

cvx_end

x_pm=x./x_env; % Phase modulation
end
