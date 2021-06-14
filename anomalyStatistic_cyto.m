close all
clear all
load('cellA.mat');
load('celldata1.mat');

x=celldata1;
N=8;

A=cellA;

%[v,d]=eig(full(A));



x=x(randperm(800,N),1:11)';
%x=x-mean(x(:,1:N/2),2);
%x=normr(x);
%x=x(:,randperm(N));
x0=x;
P=[7];
e=zeros(11,N);

ind=zeros(11,N);
ind(P,1:N)=.5*ones(length(P),N);
fas=zeros(11,N);
fas0=zeros(11,N);
for nn=1:N/2
    per=randperm(11,1);
    per=7;
x(per,nn+1:N)=5*x(per,nn+1:N);
end
for i=1:N
   % x(:,i)=x(:,i)/max(x(:,i));
    % x0(:,i)=x0(:,i)/max(x0(:,i));
    [~,~,f]=graphAS(x(:,i),A);
    
    [~,~,f0]=graphAS(x0(:,i),A);
    fas(:,i)=f;
    fas0(:,i)=f0;
    
    [fhilb(:,i),h0]=graphHT_filter(x(:,i),A,5);
    [fhilb0(:,i),h0]=graphHT_filter(x0(:,i),A,5);
end
fh=(imag(fas));
fh0=(imag(fas0));
thresh=mean(((fhilb0(:,1:N/2))),2);
thresh=kron(ones(1,N/2),thresh);
threshn=mean(((x(:,1:N/2))),2);
threshn=kron(ones(1,N/2),threshn);


figure, plot((vec(x(:,N/2+1:N))),'r'), hold on, plot((vec(x0(:,N/2+1:N))),'b'), legend('abnoral sig','normal Sig'), stem((vec(threshn)),'k');
%figure, plot(10*log10(abs(vec(fh(:,N/2+1:N)))),'r'), hold on,  plot(10*log10(abs(vec(fh0(:,N/2+1:N)))),'b'), stem(vec(ind),'k'); legend('abnoral sig','normal Sig');
figure, plot(((vec(fhilb(:,N/2+1:N)))),'r'), hold on,  plot(((vec(fhilb0(:,N/2+1:N)))),'b'), stem((vec(thresh)),'k');legend('abnoral sig','normal Sig');
%legend('abnoral sig','normal Sig','abnormal HT','normal HT');


