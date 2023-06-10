%% Main function for 2002 Sobie Simulation integrating with Breit 2018

close all;
clear;

% initialization
para1 = 0.4;   % lumenal initial condition
para2 = 1;   % V_SS
para3 = 1.0;   % \tau_{efflux}
para4 = 1;   % D_{RyR}
ryrnum = 3; % ryr num

y_Sobie_ss = mf_initialization_2002_b(para1);
tot_result = mf_equation_2002_b;

% == Calcium impulse formulation == 

I_DHPR_func = @(t)0.0-0.5*(t>0)*(t<0.5);

rosadoind = 1;
if rosadoind == 1
    j_c = 0.75;
    tau_rls = 2;
    I_DHPR_func = @(t)-1 * (t<=tau_rls).*j_c.*(1-t/tau_rls);
end

ryr_sobie_test = @(t,y) tot_result.test_ss_model(y,I_DHPR_func(t),para1,para2,para3,para4,ryrnum);

y0 = y_Sobie_ss';

opts = odeset('RelTol',1e-8,'AbsTol',1e-8,'MaxStep',1e-2);
[t,y] = ode23s(ryr_sobie_test,[0,100],y0,opts);
 
j_ryr_lst = 0*t;
for i=1:length(t)
    [dy,j_ryr] = tot_result.test_ss_model(y(i,:),I_DHPR_func(t(i)),para1,para2,para3,para4,ryrnum);
    j_ryr_lst(i) = j_ryr;
end

F = 96.485; 
V_ds = 1.0000e-13;
I_ryr_lst = 1e6 * j_ryr_lst * 2 * F * V_ds;

result = [t,y,I_ryr_lst];

save('meanfield_result_2011_alpha.mat','result');

% plot

figure()
set(gcf, 'Position',[200 200 1500 500]);

subplot(1,3,1)
plot(t,y(:,8),'LineWidth',2)
title('subspace concentration')
ssv = y(:,8);
xlim([0 15])
xlabel('Time (ms)','FontSize',16)
ylabel('Concentration (\mu M)','FontSize',16)
title('SS Concentration','FontSize',16)

subplot(1,3,2)
plot(t,y(:,7),'LineWidth',2)
xlim([0 15])
title('lumanal concentration')
xlabel('Time (ms)','FontSize',16)
ylabel('Concentration (\mu M)','FontSize',16)
title('Lumenal Concentration','FontSize',16)

subplot(1,3,3)
plot(t,y(:,2)+y(:,4),'LineWidth',2)
title('open probablity')
xlim([0 15])
xlabel('Time (ms)','FontSize',16)
ylabel('Probablity','FontSize',16)
title('Open Probablity','FontSize',16)





% figure()
% plot(t,I_ryr_lst)
% title('SR Release Flux')
% xlim([0 30])

tname = 'sobie-rosado-result';
save(['../',tname,'.mat'],'t','ssv')



