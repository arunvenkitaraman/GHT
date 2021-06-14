%% Graph Phase unwrap function

% Arun Venkitaraman 23-03-2017
% This code computes the unwrapped phase of th graph analytic signal using
% Algorithm 1 of paper
% Inputs: rawphas: phase angle of the GAS, A: normalized adj matix
% Output: phas: phase after graph unwrapping 

function [phas]=GraphUnwrap(rawphas,A)
%k Starting node
N=length(rawphas);
phas=zeros(N,1);
ord=zeros(N,1);
%for k=1:N
phas(1)=rawphas(1);
ord(1)=1;
p=1;


Omega=(1:N);
for i=2:N
    
    Omega=setdiff(Omega,p);
    zinga=zeros(1,N);
    zinga(Omega)=abs(A(p,Omega));
    [~,p]=max(zinga);
    phas(i)=rawphas(p);
    ord(i)=p;
end



phas=unwrap(phas);



