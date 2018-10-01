function EDO=HyphalTanks(t,y,N,Nv,A,D,kc,Kc,kp,Kp,m,w0,v,Yl,Yphi,Deltax,rhox,psi)

% t - time 
% y - variables to solve for: nutrient concentration, vesicle concentration, tip-tank reactor length

% 1 to N - Nutrients
% N+1 to 2*N - Vesicles
% 2*N +1 length

%last differential equation - derivative of the length of the tip tank
EDO(2*N+1)=Yl*(kc*y(2*N)/(Kc+y(2*N)));

% first set of EDOS are solver for the nutrient balance in each reactor 
% from 1 to N - Nutrient mass balance
for i=1:N
    if i == 1
       EDO(i) = (v/Deltax)*(w0-y(i))+(D/(Deltax^2))*(w0-y(i))-(D/(Deltax^2))*(y(i)-y(i+1))-m*rhox;
    elseif i >= 2 && i < (N - (Nv + 1))
       EDO(i) = (v/Deltax)*(y(i-1)-y(i))+(D/(Deltax^2))*(y(i-1)-y(i))-(D/(Deltax^2))*(y(i)-y(i+1))-m*rhox;
    elseif i >= (N - (Nv + 1)) && i <= (N - 1)
       EDO(i) = (v/Deltax)*(y(i-1)-y(i))+(D/(Deltax^2))*(y(i-1)-y(i))-(D/(Deltax^2))*(y(i)-y(i+1))...
           -(1/Yphi)*(kp*y(i)/(Kp+y(i)))-m*rhox;
    else % when i == N
       EDO(i) = -(y(i)/y(2*N+1))*EDO(2*N+1)+(v/y(2*N+1))*(y(i-1))+(D/(y(2*N+1)*(y(2*N+1)+Deltax)/2))*(y(i-1)-y(i))-m*rhox; 
    end 
end 


for i=(N+1):(2*N)
    if i == N+1
       EDO(i) = -(psi/Deltax)*y(i);
    elseif i >= N+2 && i < (2*N - (Nv+1))
       EDO(i) = (psi/Deltax)*(y(i - 1) - y(i));
    elseif i >= (2*N - (Nv + 1)) && i <= (2*N - 1)
       j= i-N; 
       EDO(i) = (psi/Deltax)*(y(i - 1)-y(i)) + (kp*y(j)/(Kp+y(j)));
    else % when i == 2*N
       EDO(i) = -(y(i)/y(2*N+1))*EDO(2*N+1) + (psi/y(2*N+1))*y(i-1) - (1/(A*y(2*N+1)))*(kc*y(i)/(Kc+y(i)));
    end 
end 

EDO = EDO';
