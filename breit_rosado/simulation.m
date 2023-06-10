%% Simulation of RyR modeling in paper: 
% Spine-to-Dendrite Calcium Modeling Discloses Relevance for Precise 
% Positioning of Ryanodine Receptor-Containing Spine Endoplasmic Reticulum

% Originally proposed in: 
% Ryanodine Receptor Adaptation and Ca2+-lnduced Ca2+ ReleaseDependent 
% Ca2+ Oscillations

% define constants
ka_minus = 28.8; % s-1
ka_plus = 1500; % \mu M-4 s-1
kb_minus = 385.9; % s-1
kb_plus = 1500; % \mu M-3 s-1
kc_minus = 0.1; % s-1
kc_plus = 1.75; % s-1

% parameter setting
rou_r = 3.0; % \mu m-2
I_ref = 3.0*10^-18; % mol s-1
btot = 40; % \mu M
c_e_ref = 200; % === "guess value", not specified in paper === 

% initial conditions
c_c_init = 0.05; % \mu M
c_e = 250; % \mu M, set as a constant throughout the process
b_init = btot*(1/(kb_plus*c_c_init/kb_minus+1))+2; % set as the steady state

% Use DHPR impulse, as mentioned in Sobie's paper
% 0.5
I_DHPR_func = @(t)0.0-1*(t>1)*(t<1.5);

% [c1,o2,c2,o1], define initial open/close probablity
% the initial conditions of parameter should be set as the final steady
% state;

% cof1 = ((ka_minus/ka_plus)/c_c^4);
% cof2 = (c_c^3)/(kb_minus/kb_plus);
% cof3 = (kc_plus/kc_minus);
% inito1 = 1/(cof1+cof2+cof3+1);
% initc1 = cof1*inito1;
% inito2 = cof2*inito1;
% initc2 = cof3*inito1;

% manually pick the steady state, after running the simulation for long time
inito1 = 0.003;
initc1 = 0.944;
inito2 = 0.000;
initc2 = 0.053;

% we should guarantee that c1,c2,o1,o2 all \ge 0
y0 = [initc1,inito2,initc2,inito1,c_c_init,b_init];
tspan = [0,10];
opts = odeset('RelTol',1e-8,'AbsTol',1e-8,'MaxStep',1e-2,'NonNegative',[1,2,3,4]);

ryr_test = @(t,y) simu_eq(y,I_DHPR_func(t));
[t,y] = ode78(ryr_test,tspan,y0,opts);

open_prob = y(:,2)+y(:,4);

close all
figure()
hold on
plot(t,y(:,2),'b') 
plot(t,y(:,4),'r')
plot(t,y(:,1),'g')
plot(t,y(:,3),'y')
plot(t,y(:,2)+y(:,4),'m','LineWidth',4)
legend('o_2','o_1','c_1','c_2','open probablity','FontSize',14)
title('RyR Open Probablity over time')

figure()
plot(t,y(:,5))
title('Cytosolic calcium concentration (c_c)')

figure()
area = 28/3; % to have roughly the same number of RyR as Sobie
mofterm = I_ref/c_e_ref; % subject to change
i_r = (c_e-y(:,5)) * mofterm;
j_r = area*rou_r.*open_prob.*i_r*10^17;
plot(t,j_r)
title('Calcium Flux')



