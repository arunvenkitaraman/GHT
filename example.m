%% Graph envelopes
% Arun Venkitaraman 23-03-2017
% This code computes the CGE and GAM for a given signal
% The same code can be used to consider
% 1. Minnesota graph
% 2. Jitter signal graph
% 3. 1D- signal graph
% 4. Social/Community graph (wieghted and unweighted)
% Copyright Arun Venkitaraman


%% Minnesota graph
%Q=load('minnesota.mat');
%A=Q.A;
%A=full(A(1:20:end,1:20:end));

%% Jitter signal ggraph
%A=circshift(diag(1+1.*rand(1,N)),1);
%%

%% 1-D signal ggraph
%A=circshift(eye(N),1);
%%


%% Community graph

Ncomm=5; % No of communities
Lcomm=10; % No of nodes within each community
N=Ncomm*Lcomm; % total no of nodes



    A=zeros(N);
    for pe=1:1:Ncomm
        A((pe-1)*Lcomm+(1:Lcomm),(pe-1)*Lcomm+(1:Lcomm))=(sprand(Lcomm,Lcomm,.5));
       
    end
   
    A=A+(1*full(sprand(N,N,0.005))); % adding sparse random edges
  
    
  


%A=sign(A); % if unweighted graph
A=A-diag(diag(A));









N=size(A,2);
n=round(N/2);
[~,d]=eig(A);
A=A/max(abs(diag(d)));
[v,d]=eig(A);
dl=diag(d);
gft=inv(v);
[bee,boo]=sort(abs(1-diag(d)/max(abs(diag(d)))));
vsort=v(:,boo);
dsort=diag(d);





%% Jitter signal

geo=randperm(3,2)+0/2;
f1=vsort(:,geo(1));
f2=vsort(:,round(N/2)+geo(2));
f=(N*real(f1.*f2));

%% Community graph signal
f=zeros(N,1);

for Com=1:2:2 % Only community 3 is active
    
    
    f(1+(Com-1)*Lcomm:(Com)*Lcomm)=(1+.1*randn(1,Lcomm));%.75*sin(100*(1:Lcomm)/Lcomm));
    
end


%   f1=vsort(:,geo(1));
%   f2=vsort(:,round(N/2)+geo(2));
%   f=(N*real(f1.*f2));
f=(f)/max(abs(f));
% f=wavread('kabree01.wav');


bet=0.1;
[am1,pm1]=CGE(f,A,bet);
[am2,pm2,fas2]=graphAS(f,A);


% Measuring MS_g for signal and envelopes
msg_f=MS_g(f,A,2)/norm(f,2)
msg_absf=MS_g(abs(f),A,2)/norm(abs(f),2)
msg_cvx=MS_g(am1,A,2)/norm(am1,2)
msg_gas=MS_g(am2,A,2)/norm(am2,2)

figure, plot(f,'b'), hold on,
plot(am1,'r'), hold on, plot(-am1,'r'), hold on,
plot(am2,'k'), hold on, plot(-am2,'k'), hold on,
legend('Signal','CGE','GAM');
