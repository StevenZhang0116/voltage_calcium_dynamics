% Initialize the data

function [y_Sobie_ss] = mf_initialization_2002(para)
    Po = 0.0;
    c_ca_lumen = 1000 * para;
    c_ca_ss = 0.1;
    
    btot = [24 47 1124];  
    kp = 1e-3 * [100 115 115];
    km = 1e-3 * [38 100 1000]; 
    b = btot.*(km./kp)./(km./kp + c_ca_ss);
        
    y_Sobie_ss = [Po,c_ca_lumen,c_ca_ss,b];
    % y_Sobie_ss = [Po,c_ca_lumen,c_ca_ss];
    
end