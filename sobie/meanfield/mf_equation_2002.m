function tot_result=mf_equation_2002
    tot_result.test_ss_model = @test_ss_model;
    tot_result.ryr_sobie = @ryr_sobie;
    tot_result.ca_er_sobie = @ca_er_sobie;
    tot_result.ca_ss_sobie = @ca_ss_sobie;
    tot_result.d_b_sobie = @d_b_sobie;
end

function [dy,j_ryr,d_y_ss,dirvv,kopenclose,facset]=test_ss_model(y,I_DHPR,para1,para2,para3,para4,ryrnum)
    % disp(y)
    % Unpack Variables
    Po = y(1); 
    c_ca_lumen = y(2);
    y_ss = y(3:end);
    c_ca_ss = y_ss(1);
    
    d_ryr = 4 * para4;
    [n_ryr,~,~] = indicator(ryrnum);
    j_ryr = n_ryr*Po*d_ryr*(c_ca_lumen-c_ca_ss); % Equation 3
    
    [d_Po,k_close,k_open,fac1,fac2] = ryr_sobie(Po,c_ca_ss,c_ca_lumen,ryrnum);
    d_ca_lumen = ca_er_sobie(c_ca_lumen,j_ryr,para1,para2);
    [d_y_ss,d1,d2,d3,d4] = ca_ss_sobie(y_ss,j_ryr,I_DHPR,para3);
    db_dt = d_b_sobie(y_ss);
    dy = [d_Po/2; d_ca_lumen; d_y_ss; db_dt(1); db_dt(2); db_dt(3)];
    dirvv = [d1,d2,d3,d4];
    kopenclose = [k_close,k_open];
    facset = [fac1,fac2];
end

function [d_Po,k_close,k_open,f1,f2] = ryr_sobie(Po,c_ca_ss,c_ca_lumen,ryrnum)
    Km_r_max = 17.14;
    alpha_r = 6.86e-3; 
    kr_plus_max = 30;
    hill = 4;
    [n_ryr,~,~] = indicator(ryrnum);
    nopen = Po*n_ryr;
    kcoop = 1;
    
    Km_r = Km_r_max-alpha_r*c_ca_lumen;
    kr_minus = 0.48; 

    % Equation 11
    kr_plus = kr_plus_max*c_ca_ss^hill/(c_ca_ss^hill+Km_r^hill); 
    f1 = c_ca_ss^hill;
    f2 = Km_r^hill;
    CF_open = 1+(nopen+1)/n_ryr; % Equation 13
    CF_closed = kcoop*(1+((n_ryr-nopen)+1)/n_ryr); % Equation 14
    k_close = kr_minus*CF_closed;
    k_open = kr_plus*CF_open;
    d_Po = (1-Po)*k_open-Po*k_close; 

end

function d_ca_lumen = ca_er_sobie(c_ca_lumen,j_ryr,para1,para2)
    CaNSR = 1000 * para1;
    tau_refill = 10;
    CSQ = 1e4;  
    KCSQ = 800;
    V_ds = 1.0000e-13 * para2;
    V_JSR = 1.0000e-11;
    
    J_refill = (CaNSR-c_ca_lumen)/tau_refill; % Equation 8
    B_JSR = (1+CSQ*KCSQ/(KCSQ+c_ca_lumen)^2)^-1; % Equation 9
    d_ca_lumen = B_JSR*(J_refill-j_ryr*V_ds/V_JSR) ; % Equation 7
end

function [d_y_ss,j_dhpr,j_efflux,j_buff,j_ryr] = ca_ss_sobie(y_ss,j_ryr,...
    I_DHPR,para3)
    F = 96485 * 1e3; 
    V_ds = 1.0000e-13;
    Camyo = 0.1;
    tau_efflux = 7e-4 * para3;
    btot = [24 47 1124]';  
    kp = 1e-3 * [100 115 115]';
    km = 1e-3 * [38 100 1000]'; 
    
    c_ca_ss = y_ss(1);
    b1 = btot.*(km./kp)./(km./kp + c_ca_ss);
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


