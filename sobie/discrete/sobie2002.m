% Regenerate the experimental simulation mentioned in Sobie, 2002. 
% Stochastic Simulation; Monte-carlo method

%% Basic Parameters
dt = 1e-5;                        % ms, adjust so that maximum transition probability = 0.01
dt_record = 0.1;
trials = 1;                       % Number of identical simulations to run
t_end = 100;
iterations = round(t_end / dt);
outputs = round(t_end / dt_record) + 1;
plottime = 0 : dt_record : (outputs - 1) * dt_record;

%% Physical Constants
F = 96.485;                       % Faraday's constant, C/mmol
% R = 8.314;                        % J mol-1 K-1
% T = 310;                          % K
% RTF = R * T / F; 

V_ds = 1.0000e-13;
V_JSR = 1.0000e-11;

tau_efflux = 7e-4;                % ms
tau_refill = 10;                  % ms

D_ryr = 4;                        % ms-1

%% RyR Parameters
global kr_minus kr_plus_max
kr_minus = 0.48;                  % ms-1 (mean open time = 2 ms) 
kr_plus_max = 30;                 % ms-1
Km_r_max = 17.14;                 % uM
alpha_r = 6.86e-3; 
hill = 4;                         % exponent

N_RyR_lst = [28];

% nstates = N_RyR + 1;
kcoop = 1;

%% Buffering parameters

% Subspace buffers, calmodulin, SR sites, SL sites
btot = [24 47 1124];             % uM
% n_buffers = length(btot);

% Factor of 1e-3 to convert from s-1 to ms-1
kp = 1e-3 * [100 115 115];       % uM^-1 ms^-1
km = 1e-3 * [38 100 1000];       % ms^-1

% JSR buffer CSQ
CSQ = 1e4;                       % uM
KCSQ = 800;                      % uM

%% Fixed ionic concentrations
Camyo = 0.01;
CaNSR = 1000;

%% Initialize arrays to hold results
Cads_all = zeros(outputs, trials);
CaJSR_all = zeros(outputs, trials);

Irel_all = zeros(outputs, trials);
Nopen_all = zeros(outputs, trials);

pmax = zeros(trials, 1); % max transition probability for adjusting dt
total_endtime = [];
total_fullopentime = [];

ratio1 = zeros(outputs, trials);

%% For loop
for i = 1 : length(N_RyR_lst)
    N_RyR = N_RyR_lst(i);
    disp(N_RyR)
    endtime_lst = [];
    fulltime_lst = [];
    
    for x = 1 : trials
        disp(x);
        % Initial conditions
        Cads = 0.1;
        CaJSR = 1000;
        nopen = 0;
        b = btot.*(km./kp)./(km./kp + Cads);
        ind = 1 ;
        count1 = 0; count2 = 0;

        for ii = 1: iterations
            time = dt * (ii - 1);
            if (time > 5 && time < 5.5)
                J_d = 1e-6 * 0.5 / (2 * V_ds * F); % Equation 2
            else
                J_d = 0 ;
            end
            J_ryr = nopen * D_ryr * (CaJSR - Cads); % Equation 3
            I_ryr = 1e6 * J_ryr * 2 * F * V_ds;
            J_efflux = (Camyo - Cads) / tau_efflux; % Equation 6
            J_refill = (CaNSR - CaJSR) / tau_refill; % Equation 8

            db_dt = -kp .* b * Cads + km .* (btot - b);
            J_buff = sum(db_dt); % Equation 5

            B_JSR = (1 + CSQ * KCSQ / (KCSQ + CaJSR)^2)^-1; % Equation 9
    
            % Write arrays after fluxes calculated, before integration and state
            if (mod(time, dt_record) == 0)
                Cads_all(ind, x) = mean(Cads);
                CaJSR_all(ind, x) = mean(CaJSR);
                Irel_all(ind, x) = I_ryr;
                Nopen_all(ind, x) = nopen;
                ind = ind + 1;
            end

            Km_r = Km_r_max - alpha_r * CaJSR;
            kr_plus = kr_plus_max * Cads^hill / (Cads^hill + Km_r^hill);
            CF_open = 1 + (nopen + 1) / N_RyR;
            CF_closed = kcoop * (1 + ((N_RyR - nopen) + 1) / N_RyR);

            % Change of amounts of open/closed channels
            pincrease = dt * (N_RyR - nopen) * kr_plus * CF_open;
            pdecrease = dt * nopen * kr_minus * CF_closed;

            if (pincrease > pmax) 
                pmax = pincrease;
            end
            if (pdecrease > pmax)
                pmax = pdecrease;
            end

            if (rand < pincrease)
                nopen = nopen + 1;
            end
            if (rand < pdecrease)
                nopen = nopen - 1;
            end


            dCads_dt = J_efflux + J_d + J_ryr + J_buff; % Equation 1
            dCaJSR_dt = B_JSR *(J_refill - J_ryr*V_ds/V_JSR) ; % Equation 7

            ratio1(index, iii) = (J_ryr * V_ds / V_JSR)/J_refill;
            

            Cads = Cads + dt*dCads_dt;
            CaJSR = CaJSR + dt*dCaJSR_dt;
            b = b + dt*db_dt;
            
            if (nopen == 0 && time > 6 && count1 == 0)
                endtime_lst(x) = time;
                count1 = 1;
            end
            
            if (nopen < 0.9*N_RyR && time > 6 && count2 == 0)
                fulltime_lst(x) = time;
                count2 = 1;
            end
        end

        % Write values at last time point t_end after loop finished
        Cads_all(end, x) = Cads;
        CaJSR_all(end, x) = CaJSR;

        Irel_all(end, x) = 1e6*J_ryr*2*F*V_ds;
        Nopen_all(end, x) = nopen;
    end
    
    avg_endtime = sum(endtime_lst) / length(endtime_lst);
    avg_fullopentime = sum(fulltime_lst) / length(fulltime_lst);
    total_endtime(i) = avg_endtime;
    total_fullopentime(i) = avg_fullopentime;
    
end


%% Average trials
Irel_all = mean(Irel_all, 2);
Cads_all = mean(Cads_all, 2); 
CaJSR_all = mean(CaJSR_all, 2);
Nopen_all_2002 = mean(Nopen_all, 2);
scale_all = rescale(Nopen_all);

%% Plot figures
% figure
% plot(plottime, Nopen_all_2002)
% title("Number of Open Channels")
% xlabel("Time")
% ylabel("Amount of Open Channels")

figure
subplot(2,2,1);
plot(plottime,Irel_all);
title("SR Ca^{2+} Release Flux");
subplot(2,2,2);
plot(plottime,Cads_all);
title("[Ca2+]_{SS}");
subplot(2,2,3);
plot(plottime,CaJSR_all);
title("[Ca2+]_{lumen}");
subplot(2,2,4);
plot(plottime,Nopen_all_2002);
title("RyR Open Probablity");

%% Save data
save('2002_plottime', 'plottime')
save('2002_numopen', 'Nopen_all_2002')
save('2002_releaseflux', 'Irel_all')
save('2002_jsr', 'CaJSR_all')

