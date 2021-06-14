%% Anomaly analysis with S_a
% Arun Venkitaraman 23-03-2017
% This code computes the S_a metric and plots them for realizations under two 
% hypothesis for anomalies
% Copyright Arun Venkitaraman
clear all
close all

Q=load('minnesota.mat');
A=full(Q.A);

N=size(A,2);
n=round(N/2);
[~,d]=eig(A);
A=A/max(abs(diag(d)));
[v,d]=eig(A);
dl=diag(d);
gft=inv(v);
[bee,boo]=sort(abs(1-diag(d)/max(abs(diag(d)))),'ascend');
vsort=v(:,boo);
dsort=diag(d);
vsort=inv(vsort);

%%
n=N;



close all
nsamp=10;
onof=[zeros(nsamp,1);ones(nsamp,1)]; % First nsamp cases belong to H1 and rest nsamp to H1
clear GFSS_sig;
clear GFSS_env;
for r=1:2*nsamp
x=zeros(n,1);
ncon=5; % no of connected points active;
for i=1:5
   
    ck=randperm(n,1);
    [y,ind]=sort(A(ck,:),'descend'); 
     
    x(ck)=5;
    
end
e=mvnrvnv(n,1);
e=pinv(.1*eye(n)+


f=onof(r)*x+e; % Observed graph signal

bet=.1;
%f=f/max(f(:));
[xe,pm1]=CGE(f,A,bet);

%figure, plot(f), hold on, plot(xe,'r'), plot(-xe,'r');


GFSS_sig(:,r)=(abs(vsort*f)).^2;

GFSS_env(:,r)=(abs(vsort*xe)).^2;

GFSS_env_sa(r)=sum(GFSS_env(1:260,r));
GFSS_sig_sa(r)=sum(GFSS_sig(1:260,r));
end
figure,plot(real(GFSS_sig_sa)),
figure, plot(real(GFSS_env_sa),'r');


