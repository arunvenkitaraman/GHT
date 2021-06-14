%% Graph AS and GHT construction
% Arun Venkitaraman 23-03-2017
% This code computes the GAS and GHT of a signal
% Inputs: f: signal, Anorm: normalized adj matix
% Output: am: graph envelope, pm: graph phase modulation, and f_as: the graph analytic signal

function [am,pm,f_gas]=graphAS(f,Anorm)


N=length(f);
[U,E]=eig(Anorm);

[bee,boo]=sort(abs(1-diag(E)/max(abs(diag(E)))));
Usort=U(:,boo);


Esort=diag(E);
Esort=Esort(boo);

F=pinv(U)*f; % FT of signal f
Fsort=F(boo); % Sorting GFT coeeficients


%% GAS spectrum creation
F_gas=0*Fsort;
for i=1:N
    
    F_gas(i)=2*Fsort(i);
        
    if (angle(Esort(i))<0&&~isequal(angle(Esort(i)),-pi)&&~isequal(angle(Esort(i)),pi))
        F_gas(i)=0;
    end
    
    if(isreal(Esort(i)))
    F_gas(i)=Fsort(i);
    end
    
end
F_gas=F_gas';


% Inverse GFT to obtain the vertex domain signal
f_gas=Usort*F_gas';
f_rec=real(f_gas); % Signal
f_hilb=imag(f_gas); % graph HT
am=abs(f_gas); % Amplitude modulation
pm=f./am; % Phase modulation
