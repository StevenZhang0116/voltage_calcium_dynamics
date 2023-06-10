% Initialize the data, Sobie 2011 + Breit

function [y_Sobie_ss] = mf_initialization_2002_b(para)
    % [c1,o2,c2,o1,c_c,b], all start from the steady state
    Po = [0.944,0.000,0.053,0.003,0.05,35.4910];
    
    c_ca_lumen = 1000 * para;
    c_ca_ss = 0.1;
    
    btot = [24 47 1124];  
    kp = 1e-3 * [100 115 115];
    km = 1e-3 * [38 100 1000]; 
    b = btot.*(km./kp)./(km./kp+c_ca_ss);
    
    % size of 11 (6+1+1+3)
    y_Sobie_ss = [Po,c_ca_lumen,c_ca_ss,b];
end