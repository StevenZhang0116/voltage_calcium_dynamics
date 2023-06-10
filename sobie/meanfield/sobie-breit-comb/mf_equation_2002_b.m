% combine with Sobie's formulation

function tot_result=mf_equation_2002_b
    tot_result.test_ss_model = @test_ss_model;
    tot_result.ryr_sobie = @ryr_sobie;
    tot_result.ca_er_sobie = @ca_er_sobie;
    tot_result.ca_ss_sobie = @ca_ss_sobie;
    tot_result.d_b_sobie = @d_b_sobie;
end

function [dy,j_ryr]=test_ss_model(y,I_DHPR,para1,para2,para3,para4,ryrnum)
    Po = y(1:6); 
    c_ca_lumen = y(7);
    y_ss = y(8:end);
    c_ca_ss = y_ss(1);

    rou_r = 5.0; % \mu m-2
    I_ref = 3.5*10^-18; % mol s-1
    c_e_ref = 250;
    
    I_r = I_ref/c_e_ref * (c_ca_lumen-c_ca_ss);
    % the parameter n_ryr is not used in the Sobie-Briet version and the
    % density function is used instead, we add a area para here to
    % match the units
    j_ryr = rou_r * (Po(2)+Po(4)) * I_r * 1e19 * 7;
    
    dydt = ryr_sobie(Po,c_ca_ss);
    d_ca_lumen = ca_er_sobie(c_ca_lumen,j_ryr,para1);
    d_y_ss = ca_ss_sobie(y_ss,j_ryr,I_DHPR,para3);
    db_dt = d_b_sobie(y_ss);
    dy = [dydt; d_ca_lumen; d_y_ss; db_dt];
end

% substitute to Breit's paper
function dydt = ryr_sobie(y,c_ca_ss)
    ka_minus = 28.8; % s-1
    ka_plus = 1500; % \mu M-4 s-1
    kb_minus = 385.9; % s-1
    kb_plus = 1500; % \mu M-3 s-1
    kc_minus = 0.1; % s-1
    kc_plus = 1.75; % s-1
    
    btot = 40; % \mu M
    kapb_plus = 27; % \mu M-1 s-1
    kapb_minus = 19; % s-1
    
    % solve eq
    dydt = zeros(4,1);
    % c_1 Eq 15
    dydt(1) = ka_minus*y(4)-ka_plus*c_ca_ss^4*y(1);
    % o_2 Eq 16
    dydt(2) = kb_plus*c_ca_ss^3*y(4)-kb_minus*y(2);
    % c_2 Eq 17
    dydt(3) = kc_plus*y(4)-kc_minus*y(3);
    % o_1 Eq 14
    dydt(4) = -(dydt(1)+dydt(2)+dydt(3));
    % c_c Eq 3, Laplacian=0
    dydt(5) = (kapb_minus*(btot-y(6))-kapb_plus*y(6)*y(5));
    % b Eq 4, Laplacian=0
    dydt(6) = (kapb_minus*(btot-y(6))-kapb_plus*y(6)*y(5));
end

function d_ca_lumen = ca_er_sobie(c_ca_lumen,j_ryr,para1)
    CaNSR = 1000 * para1;
    tau_refill = 10;
    CSQ = 1e4;  
    KCSQ = 800;
    V_ds = 1.0000e-13;
    V_JSR = 1.0000e-11;
    
    J_refill = (CaNSR-c_ca_lumen)/tau_refill; % Equation 8
    B_JSR = (1+CSQ*KCSQ/(KCSQ+c_ca_lumen)^2)^-1; % Equation 9
    d_ca_lumen = B_JSR*(J_refill-j_ryr*V_ds/V_JSR) ; % Equation 7
end

function d_y_ss = ca_ss_sobie(y_ss,j_ryr,I_DHPR,para3)
    F = 96485 * 1e3; 
    V_ds = 1.0000e-13;
    Camyo = 0.1;
    tau_efflux = 7e-4 * para3;
    btot = [24 47 1124]';  
    kp = 1e-3 * [100 115 115]';
    km = 1e-3 * [38 100 1000]'; 
    
    c_ca_ss = y_ss(1);
    b2 = y_ss(2:end);
    b2 = b2(:);

    j_dhpr = -I_DHPR/(2*F*V_ds); % Equation 2
    j_efflux = (Camyo - c_ca_ss) / tau_efflux; % Equation 6
    dbdt = -kp .* b2 * c_ca_ss + km .* (btot - b2);
    j_buff = sum(dbdt); % Equation 5
    
    d_y_ss = j_dhpr+j_efflux+j_buff+j_ryr; % Equation 1
end

function db_dt = d_b_sobie(y_ss)
    c_ca_ss = y_ss(1);
    b = y_ss(2:end);
    b = b(:);

    btot = [24 47 1124]';  
    kp = 1e-3 * [100 115 115]';
    km = 1e-3 * [38 100 1000]'; 
    
    db_dt = -kp .* b * c_ca_ss + km .* (btot - b);
end


