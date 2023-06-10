% Monte-Carlo model of Ca spark from single cluster or RyRs in Sobie, 2011.
% Recovery of calcium spark over time

%% Basic Parameters Setting
dt = 1e-5;       
dt_record = 0.1;

% Open single RyR at interval, 
% then continue to run for 'timeafter' milliseconds
interval = 5; % all second sparks terminated after a fixed duration
timeafter = 95; % time constant of recovery
trials = 1; % number of independent simulations

F = 96.485;  
V_ds = 1.0000e-12;
V_JSR = 1.6000e-12;
tau_efflux = 1.78e-3;  
tau_refill = 6.5;
EJequiv = 0.1;
D_ryr = 2.2e-12; 

%% RyR Gating Parameters
kr_minus = 0.48;               
kr_plus_max = 30;             
Km_r_max = 19.87;             
alpha_r = 1.0e-3;             
hill = 4;      
N_RyR = 28;
kcoup = exp(2 * EJequiv / (N_RyR - 1)); 


%% Buffering parameters
bt = [24 47 900];          
n_buffers = length(bt);
kp = 1e-3 * [100 115 115];     
km = 1e-3 * [38 100 1000];     

CSQ = 3e4; 
KCSQ = 630; 
Camyo = 0.1;
CaNSR = 1000;
% filename = 'result.mat';
t_end = interval + timeafter;

iterations = round(t_end / dt);
outputs = round(t_end / dt_record) + 1;

%% Vector initialization
plottime = 0 : dt_record : (outputs - 1) * dt_record;

Cads_all = zeros(outputs, trials);
CaJSR_all = zeros(outputs, trials);
Irel_all = zeros(outputs, trials);
Nopen_all = zeros(outputs, trials);
prob = zeros(outputs, trials);

total_endtime = [];
total_fullopentime = [];

ratio1 = zeros(outputs, trials);

%% For loop
% for i = 1 : length(N_RyR_lst)
%     N_RyR = N_RyR_lst(i);
    disp(N_RyR)
    endtime_lst = [];
    fulltime_lst = [];
    
    for iii= 1 : trials
        disp(iii);
        %% Initial conditions
        Cads = Camyo;
        CaJSR = CaNSR;
        nopen = 0;
        b = bt .* (km ./ kp) ./ (km ./ kp + Cads);
        index = 1;
        tlast = -dt; % time tracker
        J_d = 0;
        neverspark = 1;
        count1 = 0; count2 = 0;

        %% Start Iterations
        for ii = 1 : iterations
          time = dt * (ii - 1);
          if (time >= interval && tlast < interval)
            nopen = nopen + 5;
          end

          if (time >= interval + 10 && nopen < 1 && neverspark)
            break;
          end

          nclosed = N_RyR - nopen;

          J_ryr = nopen * D_ryr * (CaJSR - Cads) / V_ds; 
          I_ryr = 1e6 * J_ryr * 2 * F * V_ds;  
          J_efflux = (Camyo - Cads) / tau_efflux;
          J_refill = (CaNSR - CaJSR) / tau_refill;
          
          db_dt = -kp .* b * Cads + km .* (bt - b);
          J_buff = sum(db_dt);
          B_JSR = (1 + CSQ * KCSQ / (KCSQ + CaJSR) ^ 2) ^ -1;

          % Write arrays after fluxes calculated, before integration and state
          if (mod(time,dt_record) == 0)
            Cads_all(index, iii) = mean(Cads);
            CaJSR_all(index, iii) = mean(CaJSR);
            Irel_all(index, iii) = I_ryr;
            Nopen_all(index, iii) = nopen;
            index = index + 1 ;
          end

          Km_r = Km_r_max - alpha_r * CaJSR;
          kr_plus = kr_plus_max * Cads ^ hill / (Cads ^ hill + Km_r ^ hill);

          % 2011 Method
          pincrease = dt * nclosed * kr_plus * kcoup ^ (2 * nopen + 1 - N_RyR);
          pdecrease = dt * nopen * kr_minus * kcoup ^ (2 * nclosed + 1 - N_RyR);

          if (rand < pincrease)
            nopen = nopen + 1;
          end
          if (rand < pdecrease)
            nopen = nopen - 1;
          end

          if (nopen >= 10) 
            neverspark = 0;
          end

          
          if (nopen == 0 && time > 6 && count1 == 0)
              endtime_lst(iii) = time;
              count1 = 1;
          end

          if (nopen < 0.9*N_RyR && time > 6 && count2 == 0)
              fulltime_lst(iii) = time;
              count2 = 1;
          end

          % Update parameters
          dCads_dt = J_efflux + J_d + J_ryr + J_buff;
          dCaJSR_dt = B_JSR * (J_refill - J_ryr * V_ds / V_JSR); 

          ratio1(index, iii) = (J_ryr * V_ds / V_JSR)/J_refill;

          Cads = Cads + dt * dCads_dt;
          CaJSR = CaJSR + dt * dCaJSR_dt;
          b = b + dt * db_dt;

          % Update time tracker
          tlast = time;
        end

        % Write values at last time point t_end after loop finished
        Cads_all(end, iii) = Cads;
        CaJSR_all(end, iii) = CaJSR;
        Irel_all(end, iii) = 1e6 * J_ryr * 2 * F * V_ds;
        Nopen_all(end, iii) = nopen;
    end
    
    avg_endtime = sum(endtime_lst) / length(endtime_lst);
    avg_fullopentime = sum(fulltime_lst) / length(fulltime_lst);
%     total_endtime(i) = avg_endtime;
%     total_fullopentime(i) = avg_fullopentime;
% end

%% Graph Drawing
% save(filename, 'plottime', 'Cads_all', 'CaJSR_all', 'Irel_all', 'Nopen_all');

Irel_all = mean(Irel_all, 2);
Cads_all = mean(Cads_all, 2);
CaJSR_all = mean(CaJSR_all, 2);

raw_nopenall = Nopen_all;
Nopen_all_2011 = mean(Nopen_all, 2);

% raw_prob = prob;
% prob = mean(prob, 2);
% 
% figure
% plot(plottime, prob)
% xlim([0 40])

figure
plot(plottime, Nopen_all_2011)
% title("Number of Open Channels")
% xlabel("Time")
% xlim([0 40])
% ylabel("Amount of Open Channels")
% 
save("2011_numopen", 'Nopen_all_2011')

sparkdices = find(max(Nopen_all_2011) > 5);

figure
plot(plottime,CaJSR_all(:,sparkdices))
hold on 
CaJSR_avg = mean(CaJSR_all(:, sparkdices), 2);
plot(plottime, CaJSR_avg, 'b', 'LineWidth', 2.25)
print -depsc freeCaJSR

CaJSRtot_all = CaJSR_avg + ...
  CaJSR_avg*CSQ./(KCSQ + CaJSR_avg) ;
CaJSRtot_all = CaJSRtot_all/max(CaJSRtot_all);
figure
plot(plottime,CaJSRtot_all)
hold on
plot(plottime, 1 - exp(-plottime/90),'r','LineWidth',2.5)
xlabel("Time after Ca^{2+} spark initiation (ms)")
ylabel("Normalized total JSR [Ca^{2+}]")
print -depsc totalCaJSR

save('2011_jsr', 'CaJSRtot_all', 'CaJSR_avg');

figure
plot(plottime, Irel_all)
% xlim([0 40])
save('2011_releaseflux', 'Irel_all');

