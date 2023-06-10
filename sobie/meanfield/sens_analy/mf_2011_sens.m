%% Main function for 2011 Sobie Simulation
%% Including several calculations of sensitivity analysis

clear;

tic
breakval = 199;

variation1 = linspace(0.1,10,breakval);
% variation1 = [0.1,1,5,10];
% variation1 = [1];

% variation1 = linspace(0.01,0.095,18);
% variation2 = linspace(0.1,0.95,18);
% variation3 = linspace(1,9.5,18);
% varitaion4 = linspace(10,100,19);
% variation5 = linspace(200,1000,9);
% variation1 = [variation1,variation2,variation3,varitaion4,variation5];
% variation2 = linspace(1,100,breakval*10);
% variation = [variation1,variation2]; % coefficients use to modify parameter
% variation1 = [1];
% 
% 
variation = variation1;

% paraname = 'paras2.mat'; 
% para_generator(paraname,0)
% 
% variation = load('paras2.mat').all_parameters;

all_result = zeros(1,7);
% all_result = zeros(1,6);
para = constants.para;

stable_lst = zeros(1,11);

folder = "./data"+constants.para+"/";
name_prefix = folder+"result";

% figure()
% hold on;

lastresult = 0;

% the order of parameters: 
% alpha, dryr, nryr, csq, efflux, refill, vjsr, vss, bsl, kcsq, kmax
for jj = 1:length(variation)
%     disp(jj)
%     coeff = variation(jj,:);
    coeff = variation(jj);
    disp(coeff)
    all_result(jj,1) = coeff;

    %% demo test
%     coeff = [1.0000    1.0500    0.9500    0.9000    1.0000    1.0500    1.0500    1.0000    1.0000    0.9500    1.0000];
%     disp(coeff)

    check = ones(11);
    check(10) = coeff;
    coeff = check;

    %% initialization
    mf_initialization_2011(coeff(9)) % change to coeff if modify BSL total in sensitivity analysis test

    tot_result = mf_equation_2011;
    
    if constants.ind == 0
        I_DHPR_func = @(t)0.0-0.5*(t>5)*(t<5.5);
    else
        I_DHPR_func = @(t)0.0;
    end
    
    bump_func = @(t)0; %% save here for reference, might be useful
    
    ryr_sobie_test = @(t,y) tot_result.test_ss_model(y,I_DHPR_func(t),bump_func(t),coeff);
    
    load('y_Sobie_ss_2011.mat')
    y0 = y_Sobie_ss';
    
    opts = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',1e-1);
    [t,y] = ode45(ryr_sobie_test,[0,100],y0,opts);
     
    j_ryr_lst = 0*t;
    for i=1:length(t)
        [dy,j_ryr] = tot_result.test_ss_model(y(i,:),I_DHPR_func(t(i)),bump_func(t(i)),coeff);
        j_ryr_lst(i) = j_ryr;
    end
    
    % figure();
    % plot(t,y(:,1));
    
    F = 96.485; 
    V_ds = 1.0000e-12;
    I_ryr_lst = 1e6 * j_ryr_lst * 2 * F * V_ds;
    % figure()
    % plot(t,I_ryr_lst);
    
    if constants.ind == 1
        t = t + 5; %% push back
        k = length(t)/20;
        t_init = linspace(0,5,k)';
        result = [t,y,I_ryr_lst];
        zeromat = zeros(length(t_init),1);
        jsrarray = ones(length(t_init),1) * 1000;
        initmat = [t_init,zeromat,jsrarray,zeromat,zeromat,zeromat,zeromat,zeromat];
        check = [initmat;result];
        
        cutoff = 0;
        for i=1:length(check(:,1))
            if check(i,1) > 100
                cutoff = i;
                break;
            end
        end
        
        check = check(1:cutoff, :);
        result = check;
    else
        result = [t,y,I_ryr_lst];
    
    end
    
%     plot(result(:,1),result(:,2),'LineWidth',1.5);
%     plot(result(:,1),result(:,3),'LineWidth',1.5);
%     plot(result(:,1),result(:,8),'LineWidth',1.5);
%     scatter(log(coeff), max(result(:,8)),'b');

    %% save data
    if(~exist(folder,'dir'))
        mkdir(folder);
    end

    tempname = name_prefix+int2str(jj);
    save(tempname,"result");

    name = "meanfield_result_2011_" + para;
    save(name,'result');
    
%     figure()
%     plot(result(:,1),result(:,8))

    %% run raw sensitivity analysis
%     [result1, sens] = raw_sens_analysis(name, coeff);
%     all_result(jj,2:end) = result1;
%     if (sens ~= 0)
%         stable_lst(end+1,:) = sens; 
%         all_result(end+1,:) = result1;
%         disp(sens)
%     end
end

% save('result2.mat','stable_lst','all_result');

% xlabel('k')
% ylabel('max val')
% title('Global Maximum of SR Release Flux under different k')

% xlabel('t') 
% title('Normalized open probablity under different k')
% title('Lumen concentration under different k')
% title('SR Release flux under different k')
% legend('k=0.1','k=1','k=5','k=10');
% hold off

%% rela sens = 1/4(x(k-2)+x(k-1)+x(k+1)+x(k+2))/k ~ central difference
for jj = 3:length(variation)-2
    disp(jj);
    coeff = variation(jj);
    before_data1 = load(name_prefix+int2str(jj-1)).result;
    before_data2 = load(name_prefix+int2str(jj-2)).result;
    curr_data = load(name_prefix+int2str(jj)).result;
    next_data1 = load(name_prefix+int2str(jj+1)).result;
    next_data2 = load(name_prefix+int2str(jj+2)).result;
    avg_result = 1/4*(rela_analysis(curr_data,before_data1)+rela_analysis(curr_data,before_data2)+ ...
        rela_analysis(curr_data,next_data1)+rela_analysis(curr_data,next_data2));
    all_result(jj,2:end) = avg_result;
end

dataname = "vari-"+para;
save(dataname,"all_result");
tlst = raw_sens_graph(dataname);
disp(tlst);
toc
