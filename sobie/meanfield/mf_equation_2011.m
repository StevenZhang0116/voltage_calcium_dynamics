%% Functions

function tot_result=mf_equation_2011
    tot_result.test_ss_model = @test_ss_model;
    tot_result.ryr_sobie = @ryr_sobie;
    tot_result.ca_er_sobie = @ca_er_sobie;
    tot_result.ca_ss_sobie = @ca_ss_sobie;
    tot_result.d_b_sobie = @d_b_sobie;
end

% the order of parameters: 
% alpha, dryr, nryr, csq, efflux, refill, vjsr, vss, bsl, kcsq, kmax
function [dy,j_ryr]=test_ss_model(y,I_DHPR,para2,para3)
    % disp(y)
    % Unpack Variables
    Po = y(1); 
    c_ca_lumen = y(2);
    y_ss = y(3:end);
    c_ca_ss = y_ss(1);
    
    d_ryr = 2.2;
    n_ryr = constants.n_ryr;
    j_ryr = n_ryr*Po*d_ryr*(c_ca_lumen-c_ca_ss); 
    
    d_Po = ryr_sobie(Po,c_ca_ss,c_ca_lumen);
    d_ca_lumen = ca_er_sobie(c_ca_lumen,j_ryr,para2);
    d_y_ss = ca_ss_sobie(y_ss,j_ryr,I_DHPR,para3);
    db_dt = d_b_sobie(y_ss);
    dy = [d_Po; d_ca_lumen; d_y_ss; db_dt(1); db_dt(2); db_dt(3)];
end

function d_Po = ryr_sobie(Po,c_ca_ss,c_ca_lumen)
    Km_r_max = 19.87;
    alpha_r = 1.0e-3; 
    kr_plus_max = 30;
    hill = 4;
    n_ryr = constants.n_ryr;
    EJequiv = 0.1;

    nopen = Po*n_ryr;
    nopen = nopen;
    nclosed = n_ryr-nopen;

    kcoup = exp(2*EJequiv/(n_ryr - 1));
    kr_minus = 0.48;

    Km_r = Km_r_max-alpha_r*c_ca_lumen;
    kr_plus = kr_plus_max*c_ca_ss^hill/(c_ca_ss^hill+Km_r^hill); 
    
    d_Po = (1-Po)*kr_plus*kcoup^(2*nopen+1-n_ryr)-Po*kr_minus...
        *kcoup^(2*nclosed+1-n_ryr);
end

function d_ca_lumen = ca_er_sobie(c_ca_lumen,j_ryr,para2)
    CaNSR = 1000;
    tau_refill = 6.5;
    CSQ = 3e4;  
    KCSQ = 630;
    V_ds = 1.0000e-12 * para2;
    V_JSR = 1.6000e-12;
    
    J_refill = (CaNSR-c_ca_lumen)/tau_refill; 
    B_JSR = (1+CSQ*KCSQ/(KCSQ+c_ca_lumen)^2)^-1; 
    d_ca_lumen = B_JSR*(J_refill-j_ryr*V_ds/V_JSR); 
end

function d_y_ss = ca_ss_sobie(y_ss,j_ryr,I_DHPR,para3)
    F = 96485 * 1e3; 
    
    V_ds = 1.0000e-12;
    Camyo = 0.1;
    tau_efflux = 1.78e-3*para3;
    btot = [24 47 900]';  
    kp = 1e-3 * [100 115 115]';
    km = 1e-3 * [38 100 1000]'; 
    
    c_ca_ss = y_ss(1);
    b2 = y_ss(2:end);
    b2 = b2(:);

    j_dhpr = -I_DHPR/(2*F*V_ds); 
    j_efflux = (Camyo - c_ca_ss) / tau_efflux; 
%     disp("j_efflux"+j_efflux)
    dbdt = -kp .* b2 * c_ca_ss + km .* (btot - b2);
    j_buff = sum(dbdt);
    
    d_y_ss = j_dhpr+j_efflux+j_buff+j_ryr; 
end

function db_dt = d_b_sobie(y_ss)
    c_ca_ss = y_ss(1);
    b = y_ss(2:end);
    b = b(:);

    btot = [24 47 900]';  
    kp = 1e-3 * [100 115 115]';
    km = 1e-3 * [38 100 1000]'; 

    db_dt = -kp .* b * c_ca_ss + km .* (btot - b);
end


