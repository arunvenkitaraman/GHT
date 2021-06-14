%% Mean square measure for graph signals
% Arun Venkitaraman 23-03-2017
%p:=norm space 1 or 2;
function t=MS_g(x,Anorm,p)

t=norm(x-Anorm*x,p);%./norm(x,p);