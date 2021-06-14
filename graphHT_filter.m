%% Graph Hilbert transform as a graph filter
% Arun Venkitaraman 23-03-2017
% This code computes the GHT of a signal as a graph filter through
% computing filter coefficients
% Inputs: x: signal, Anorm: normalized adj matrix, L: # filter taps
%
% Output: xhilb, GHT of x, and corresponding filter taps h0
function [fhilb,h0]=graphHT_filter(f,Anorm,L)
% L = no of filter taps

A=Anorm;
n=size(A,2);
%L=10;
[v,d]=eig(A);
b=zeros(n,1);
lam=diag(d);
A=A/max(abs(lam));
[v,d]=eig(A);
lam=diag(d);


van_lam=vander(lam);
van_lam=fliplr(van_lam(:,end-L+1:end));
b=-1i.*sign(angle(lam));
b(angle(lam)==0)=0;


% Least squares solution for filter taps
cvx_begin quiet
variable h0(L,1) 
minimize(norm(van_lam*h0-b,2)+.001*norm(h0,2))
cvx_end 


% generating the GHT as a filter output
fhilb=0;
temp=eye(n);
for i=1:L
fhilb=fhilb+real((h0(i)))*temp*f;
temp=temp*A;
end


