%% Main function for 2011 Sobie Simulation integrating with Breit 2018

close all;
clear;

%% initialization
mf_initialization_2011_b() 

tot_result = mf_equation_2011_b;

ind = 0;

if ind == 0
    I_DHPR_func = @(t)0.0-0.5*(t>5)*(t<5.5);
else
    I_DHPR_func = @(t)0.0;
end

ryr_sobie_test = @(t,y) tot_result.test_ss_model(y,I_DHPR_func(t));

load('y_Sobie_ss_2011_b.mat')
y0 = y_Sobie_ss';

opts = odeset('RelTol',1e-8,'AbsTol',1e-8,'MaxStep',1e-2);
[t,y] = ode23s(ryr_sobie_test,[0,100],y0,opts);
 
j_ryr_lst = 0*t;
for i=1:length(t)
    [dy,j_ryr] = tot_result.test_ss_model(y(i,:),I_DHPR_func(t(i)));
    j_ryr_lst(i) = j_ryr;
end

F = 96.485; 
V_ds = 1.0000e-12;
I_ryr_lst = 1e6 * j_ryr_lst * 2 * F * V_ds;

result = [t,y,I_ryr_lst];

save('meanfield_result_2011_alpha.mat','result');

figure()
plot(t,y(:,2)+y(:,4))
title('open probablity')
figure()
plot(t,y(:,7))
title('lumanal concentration')
figure()
plot(t,I_ryr_lst)
title('SR Release Flux')


