%% Simple simulation code for Hodgkin-Huxley Model

close all
clear 
clc

simulationTime = 100;
deltaT = 0.01;
t = 0:deltaT:simulationTime;

changeTimes = [0];
currentLevels = [3];

I(1:500) = currentLevels;
I(501:2000) = 0;
I(2001:numel(t)) = currentLevels;

% Table 3
gbar_K = 36; gbar_Na = 120; gbar_L = 0.3;
E_K = -12; E_Na = 115; E_L = 10.6;
C = 1;

V = 0;
alpha_n = 0.01*((10-V)/(exp((10-V)/10)-1)); % Equation 12
beta_n = 0.125 * exp(-V/80); % Equation 13
alpha_m = 0.1*((25-V)/(exp((25-V)/10)-1)); % Equation 20
beta_m = 4*exp(-V/18); % Equation 21
alpha_h = 0.07*exp(-V/20); % Equation 23
beta_h = 1/(exp((30-V)/10)+1); % Equation 24

n(1) = alpha_n/(alpha_n+beta_n); % Equation 9
m(1) = alpha_m/(alpha_m+beta_m); % Equation 18
h(1) = alpha_h/(alpha_h+beta_h); % Equation 18

for i=1:numel(t)-1
    alpha_n(i) = 0.01*((10-V(i))/(exp((10-V(i))/10)-1));
    beta_n(i) = 0.125 * exp(-V(i)/80);
    alpha_m(i) = 0.1*((25-V(i))/(exp((25-V(i))/10)-1));
    beta_m(i) = 4*exp(-V(i)/18);
    alpha_h(i) = 0.07*exp(-V(i)/20);
    beta_h(i) = 1/(exp((30-V(i))/10)+1);

    I_Na = (m(i)^3)*gbar_Na*h(i)*(V(i)-E_Na); % Equation 3 and 14
    I_K = (n(i)^4)*gbar_K*(V(i)-E_K); % Equation 4 and 6
    I_L = gbar_L*(V(i)-E_L); % Equation 5
    I_ion = I(i)-I_K-I_Na-I_L;

    V(i+1) = V(i)+deltaT*I_ion/C;
    n(i+1) = n(i)+deltaT*(alpha_n(i)*(1-n(i))-beta_n(i)*n(i)); % Equation 7
    m(i+1) = m(i)+deltaT*(alpha_m(i)*(1-m(i))-beta_m(i)*m(i)); % Equation 15
    h(i+1) = h(i)+deltaT*(alpha_h(i)*(1-h(i))-beta_h(i)*h(i)); % Equation 16

end

V = V-70;
plot(t,V,'LineWidth',3)
hold on
legend({'voltage'})
ylabel('Voltage(mv)')
xlabel('time(ms)')
title('Voltage over time in simulated neuron')

figure()
p1 = plot(t,gbar_K*n.^4,'LineWidth',2);
hold on
p2 = plot(t,gbar_Na*(m.^3).*h,'r','LineWidth',2);
legend([p1,p2],'Potassium Conductance','Sodium Conductance')
ylabel('Conductance')
xlabel('time(ms)')
title('Conductance for Potassium and Sodium Ions in Simulated Neruon')

figure()
hold on
plot(t,n,'LineWidth',3)
plot(t,m,'LineWidth',3)
plot(t,h,'LineWidth',3)
legend({'n','m','h'})
hold off




