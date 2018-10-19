% A Model for Growth of a Single Fungal Hypha Based on Well-Mixed Tanks
% in Series: Simulation of Nutrient and Vesicle Transport in Aerial Reproductive Hyphae 
% AUTHORS: Balmant, W. et al., PLoS One 10, e0120307 (2015).
% ________________________________________________________________________
% Implementation by: Leandro A. Scholz (leandro dot a dot scholz at gmail dot com)
% MSc. in Chemical Engineering
% Federal University of Paraná
% Date: September 28th 2018
% 
% Comment: this version solves the model for N taaks and finds the time
% elapsed for the tip tank to reach twice its length.
% ________________________________________________________________________
% KEY INPUT PARAMETERS
%
% N - Initial number of tanks
% Nv - Number of vesicle producing tanks     
N = 50;
Nv = 30; 
        
% OTHER INPUT PARAMETERS 
% 
% A      - Cross sectional area of the hypha - dm^2 
A = 1e-8; 
% A0     - Cross sectional area of the source tank - dm^2
% D      - Difusivity of nutrient inside the hypha - dm^2.h^-1
D = 2.48e-4;
% kc     - maximum rate of vesicle consumption - g-vesicles.h^-1
kc = 2e-8;
% Kc     - saturation constant for vesicle consumption - g-vesicles.dm^-3
Kc = 400;
% kp     - Maximum rate of vesicle production - g-vesicles.dm^-3.h^-1
kp = 1000; 
% Kp     - saturation constant for vesicle production - g-vesicles.dm^-3
Kp = 10;
% m      - maintenance coefficient of the hypha for nutrient - g-nutrient.g-biomass^-1.h^-1
m = 1.8e-3;
% w0     - concentration of nutrient in the source tank - g-nutrient.dm^-3 
w0 = 5; 
% v      - convective velocity inside the hypha - dm.h^-1 
v = 0.026;
% Yl     - extension of hyphal length per mass of vesicles consumed - dm.g-vesicles^-1  
Yl = 1e6; 
% Yphi   - Yield coefficient for production of vesicles from nutrient - g-vesicles.g-nutrient^-1
Yphi = 0.5; 
% Deltax - length of the side of each cubic tank - dm 
global Deltax
Deltax = 1e-4; 
% to be used in Events function 
% lambda  - Maximum possible length of the vesicle producing zone - micrometers 
lambda = Nv * Deltax;
% rhox   - biomass dry weight per volume - g-biomass.dm^-3 
rhox = 100; 
% psi    - velocity of active transport of vesicles inside the hypha - dm.h^-1 
psi = 0.05; 
%   
%_________________________________________________________________________
% VARIABLES
%       
% L - length of the tip tank
% wi - concentration of nutrient in tank i
% phi_i - concentration of vesicles in tank i 
% n - number of tanks present in the hypha at any time at time t=0 -> n=N
% The number of equations to solve depends on n 
% Number of equations = 2*n + 1

% anonymous function that is just a function of t and the variables to solve
% for, not the input parameters. 
F=@(t,y) HyphalTanks(t,y,N,Nv,A,D,kc,Kc,kp,Kp,m,w0,v,Yl,Yphi,Deltax,rhox,psi);

% initial conditions 
y0=zeros(1,2*N+1); 
y0(end)= Deltax; 

%time span 
tspan=0:0.01:10;

%set Event function in options
options = odeset('Events',@myEvent);

[t,y,te,ye,ie]=ode45(F,tspan, y0,options);
