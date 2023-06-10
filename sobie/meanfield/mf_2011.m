%% Main function for 2011 Sobie Simulation

close all;
clear;

%% initialization

param1 = 0.2:0.01:0.2;
param2 = 0.7:0.5:0.7;
result1 = [];
result2 = [];

for ii = 1:length(param1)
    disp(param1(ii))
    for jj = 1:length(param2)
        para1 = param1(ii);
        para2 = param2(jj); 
        disp(para2)

        para3 = 1; % efflux
        
        mf_initialization_2011(para1); 
    
        tot_result = mf_equation_2011;
        
        if constants.ind == 0
            I_DHPR_func = @(t)0.0-0.5*(t>5)*(t<5.5);
        else
            I_DHPR_func = @(t)0.0;
        end
        
        bump_func = @(t)0; %% save here for reference, might be useful
        
        ryr_sobie_test = @(t,y) tot_result.test_ss_model(y,I_DHPR_func(t),para2,para3);
        
        load('y_Sobie_ss_2011.mat')
        y0 = y_Sobie_ss';
        
        opts = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',1e-1);
        [t,y] = ode45(ryr_sobie_test,[0,100],y0,opts);
         
        j_ryr_lst = 0*t;
        for i=1:length(t)
            [dy,j_ryr] = tot_result.test_ss_model(y(i,:),I_DHPR_func(t(i)),para2,para3);
            j_ryr_lst(i) = j_ryr;
        end
        
        % figure();
        % plot(t,y(:,1));
        
        F = 96.485; 
        V_ds = 1.0000e-12;
        I_ryr_lst = 1e6 * j_ryr_lst * 2 * F * V_ds;
        % figure()
        % plot(t,I_ryr_lst);
        
        result = [t,y,I_ryr_lst];
        
        save('meanfield_result_2011_alpha.mat','result');
         
        prob = y(:,1);
        cal_ss = y(:,3);
    
        test
    
        result1(end+1,1) = diff1; result1(end,2) = diff11; result1(end,3) = para1; result1(end,4) = para2;
        result2(end+1,1) = diff2; result2(end,2) = diff22; result2(end,3) = para1; result2(end,4) = para2;
    end
end

% figure()
% plot3(result1(:,3),result1(:,4),result1(:,1))
% xlabel('coeff of change in [Ca^{2+}]_{lumen} initial condition')
% ylabel('coeff of change in V_{SS}')
% zlabel('diff of maximum value of open probability')
% title('change in open probability')
% 
% figure()
% plot3(result2(:,3),result2(:,4),result2(:,1))
% xlabel('coeff of change in [Ca^{2+}]_{lumen} initial condition')
% ylabel('coeff of change in V_{SS}')
% zlabel('diff of maximum value of calcium concentration in SR')
% title('change in SR concentration')

% figure()
% plot(param,result1(:,1))
% hold on
% % plot(param,result1(:,2)/1000)
% xlabel('param')
% ylabel('amount of changing')
% legend('max difference','integral difference (\times 10^3)')
% title(['max open prob change, original max value = ',num2str(max(std1))])
% hold off
% 
% figure()
% plot(param,result2(:,1))
% hold on
% % plot(param,result2(:,2)/1000)
% xlabel('param')
% ylabel('amount of changing')
% legend('max difference','integral difference (\times 10^3)')
% title(['max c_c calcium concentration change, original max value = ',num2str(max(std2))])
% hold off

save('result1.mat','result1');
save('result2.mat','result2');

%% analyze result
aaaa = load("result1.mat").result1;
bbbb = load("result2.mat").result2;
ind = aaaa(:,1) < 0.8;
tt = bbbb(ind,1);
[val,idx] = max(tt);
k=0;
for i=1:length(ind)
    if ind(i) == 1
        k = k+1;
    end
    if k==idx
        output = i;
        break
    end
end
t1 = aaaa(output,3);
t2 = aaaa(output,4);
t1
t2





