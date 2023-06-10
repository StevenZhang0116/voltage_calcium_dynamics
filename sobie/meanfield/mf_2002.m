%--------------------------------------------------------------------------
% Mean-field Representation of Sobie 2002 Model
% Main Function
% Many thanks to Guanchun Li, Courant Institute for his generous assistence 
% and inspiration
%
% Steven Zhang, Courant Institute
% Updated July 2022
%--------------------------------------------------------------------------

clear;
% clc
close all
 
param1 = 0.95:0.05:0.95; % lumenal initial condition
param2 = 1:0.5:1;     % V_SS
param3 = 0.4:0.05:0.4; % \tau_{efflux}
param4 = 1.25:0.01:1.25; % D_{RyR}
ryrnumpara = 10:4:10; % ryr num
M1 = combvec(param1,param2,param3,param4,ryrnumpara);

tt = size(M1,2);
result1 = zeros(tt,1);
result2 = zeros(tt,1);
result3 = zeros(tt,1);
result4 = zeros(tt,1);

tttt = load('./res/test.mat');

std1 = load('stdopenprob.mat').prob;
std2 = load('stdcalss.mat').calss;
std3 = load('stdcallumen.mat').stdcallumen;
stdt = load('stdtime.mat').t;

for iii = 1:tt
    disp(iii)
    MM = M1(:,iii);

    para1 = MM(1);
    para2 = MM(2);
    para3 = MM(3);
    para4 = MM(4);
    ryrnum = MM(5);

    y_Sobie_ss = mf_initialization_2002(para1);
    tot_result = mf_equation_2002;
    
    % Original DHPR
    I_DHPR_func = @(t)0.0-0.5.*(t>=5).*(t<=5.5);

    % ====================
    % DHPR formulation
    % ====================  
    
    % whether to smooth using Gaussian
    smoothind = 0;
    convind = 0;
    if smoothind == 1
        int = 0.001; % length of interval
        mm = 5.25;
        sig = 0.15; % sigma, as the standard deviation
        % Gaussian function
        adj = 5*abs(sig-0.1)*(sig>0.1);
        origtm = 5;
        tm = 5-adj;
        tp = 5.5+adj;
        ss = 4 - adj;
        ee = 6 + adj;
        xx = ss:int:ee;
    
        modfy_DHPR = @(t)((1/(sig*sqrt(2*pi)).*exp(-(t-mm).^2/...
            (2*sig^2)))/4 *-1);
    
        summ1 = [];
        for hh=1:length(xx)
            summ1(end+1) = modfy_DHPR(xx(hh));
        end
        sum(summ1) * int % integral of gaussian

        % original result
        figure()
        plot(xx,modfy_DHPR(xx))
        hold on
        plot(xx,I_DHPR_func(xx))
        hold on
        legend('gaussian','step')
        title('gaussian')
        pause
        
        if smoothind == 1 && convind ~= 1
            I_DHPR_func = modfy_DHPR;
        end
        
        if convind == 1
            res1 = I_DHPR_func(xx);
            res2 = modfy_DHPR(xx);
            cr = conv(res1,res2,'same');
            aa = sum(cr)*2;
            % convolution result
            adfac = -500/aa;
            cr = cr * adfac;
        
            sum(cr) * int % integral of convolution
        
            newxx = xx-(mm-origtm);
        
            % convolution result
            figure()
            hold on
            plot(newxx,cr)
            plot(xx,I_DHPR_func(xx))
            
            title('convolution')
            pause;
        
            fz = 0;
            sz = 0;
            for nn=1:length(newxx)
                ll = cr(nn);
                if round(ll,4) < -0.0001 
                    fz = nn;
                    break
                end
            end
        
            for nn=1:length(newxx)
                ll = cr(nn);
                if round(ll,4) > -0.0001 && newxx(nn) > mm
                    sz = nn;
                    break
                end
            end
        
            % piecewise interpolation
            pip = interp1(newxx,cr,'linear','pp');
        
            plot(newxx,ppval(pip,newxx))
            legend('convolution','step','linear interpolation')
            hold off
            pause;
        
            I_DHPR_func = @(t) (ppval(pip,t) .* (t>=newxx(1)) .* (t<=newxx(end)) + 0);
        end
    end

    % ====================
    % start calculation
    % ====================  
    close all

    ryr_sobie_test = @(t,y) tot_result.test_ss_model(y,I_DHPR_func(t),...
        para1,para2,para3,para4,ryrnum);
   
    y0 = y_Sobie_ss';
    % y0(1) = 0;
    
    % MOD by Li 6
    opts = odeset('RelTol',1e-6,'AbsTol',1e-7,'MaxStep',1e-1);
    [t,y] = ode45(ryr_sobie_test,[0,100],y0,opts);
    
    j_ryr_lst = 0*t;
    dlst = zeros(length(t),4);
    dyss_lst = 0*t;
    kcloseset = 0*t;
    kopenset = 0*t;
    backval = 0;
    f1set = 0*t; 
    f2set = 0*t;
    for ii=1:length(t)
        [dy,j_ryr,d_y_ss,dirvv,kopenclose,facset] = tot_result.test_ss_model(...
            y(ii,:),I_DHPR_func(t(ii)),para1,para2,para3,para4,ryrnum);
        j_ryr_lst(ii) = j_ryr;
        dlst(ii,1) = dirvv(1);
        dlst(ii,2) = dirvv(2);
        dlst(ii,3) = dirvv(3);
        dlst(ii,4) = dirvv(4);
        dyss_lst(ii) = d_y_ss;
        kcloseset(ii) = kopenclose(1);
        kopenset(ii) = kopenclose(2);
        f1set(ii) = facset(1);
        f2set(ii) = facset(2);
    end

%     figure()
%     for ind = 1:4
%         hold on
%         plot(t,dlst(:,ind),'LineWidth',2);
%     end
%     plot(t,dyss_lst,'LineWidth',2,'LineStyle','-.')
%     hold off
%     legend('j_{dhpr}','j_{efflux}',',j_{buff}','j_{ryr}','dy_{ss}')
%     xlim([0 20])
%     xlabel('time')


    F = 96.485; 
    V_ds = 1.0000e-13;
    I_ryr_lst = 1e6 * j_ryr_lst * 2 * F * V_ds;
    
    result = [t,y,I_ryr_lst];
    
    prob = y(:,1);
    cal_ss = y(:,3);
    cal_lumen = y(:,2);
    flux = I_ryr_lst;

    % ===================================== % 
    % calculate ratios
    % ===================================== % 
    
    % calculate the maximum difference and also the integral difference
    diff1 = max(prob);
    diff11 = sum(prob);
    diff2 = max(cal_ss);
    diff22 = sum(cal_ss);
    diff23 = sum(cal_ss) * (max(t)/length(t));

    probendtime = 0;
    
    for i=1:length(stdt)
        if stdt(i) > 6
            ind = i;
            break
        end
    end
    
    for j=ind:length(stdt)
        if std1(j) < 0.001
            stdendtime = stdt(j);
            break
        end
    end
    
    for i=1:length(t)
        if t(i) > 6
            ind2 = i;
            break
        end
    end
    
    for j=ind2:length(t)
        if prob(j) < 0.01
            probendtime = t(j);
            break
        end
    end
    
    timediff = probendtime;
%     timediff

    
    hold off

    figure()
    hold on
    plot(t,prob,'LineWidth',2)
    plot(tttt.t,tttt.prob,'-.','LineWidth',1)
    ylim([0 1])
    legend('smooth','orig')
    title('open prob comparsion')
    xlabel('time (ms)')
    ylabel('open probablity')
    hold off
    xlim([0 14])
    

    figure()
    hold on
    plot(t,cal_lumen,'LineWidth',2)
    plot(tttt.t,tttt.cal_lumen,'-.','LineWidth',1)
    legend('smooth','orig')
    title('ER calcium concentration comparsion')
    xlabel('time (ms)')
    ylabel('concentration (\mu M)')
    ylim([0 1000])
    hold off

    figure()
    hold on
    plot(t,cal_ss,'LineWidth',2)
    plot(tttt.t,tttt.cal_ss,'-.','LineWidth',1)
    legend('smooth','orig')
    title('SR calcium concentration comparsion')
    xlabel('time (ms)')
    ylabel('concentration (\mu M)')
    xlim([0 14])

    figure()
    plot(t,kcloseset)
    hold on
    plot(t,kopenset)
    hold off
    xlabel('time (ms)')
    ylabel('rate of open/close ms^{-1}')
    xlim([0 14])
    ylim([0 20])
    legend('k_{close}','k_{open}')


    % open prob diff
    result1(iii) = diff1; 
    % calcium maximum diff
    result2(iii) = diff2; 
    % complete closing time diff
    result3(iii) = timediff;

    for kk = 1:length(t)
        tt = t(kk);
        if round(tt,2) == 5.5
            backind = kk;
        end
    end
    [firstmax,t1] = max(cal_ss(1:backind));
    [secondmax,t2] = max(cal_ss(backind:end));
    t2 = backind + t2;
    vdiff1 = firstmax-secondmax;
    vdiff2 = firstmax-min(cal_ss(t1:t2));
    maxtimediff = t(t2)-t(t1);
    result4(iii) = maxtimediff;
end

% ===================================== % 
% save data
% ===================================== % 

% save('result1','result1')
% save('result2','result2')
% save('result3','result3')
% save('result4','result4')
% save('M1','M1')

if smoothind == 1 && convind == 1
    name = ['smooth-',num2str(sig)];
elseif smoothind == 1 && convind ~= 1
    name = ['gaussian-',num2str(sig)];
else
    name = ['test'];
end

save(['./res/',name,'.mat'],'t','cal_lumen','prob','I_DHPR_func','cal_ss')


% ===================================== % 
% find the optimal result
% ===================================== % 

% con = load('result2.mat').result2;
% time = load('result3.mat').result3;
% maxdiff = load('result4.mat').result4;
% M1 = load('M1.mat').M1;
% 
% small = 1e10;
% ind = 0;
% for jjjj = 1:length(con)
%     if con(jjjj) < small && maxdiff(jjjj) > 0.01 && time(jjjj) < 15
%         small = con(jjjj);
%         ind = jjjj;
%     end
% end
% 
% M1(:,ind)






