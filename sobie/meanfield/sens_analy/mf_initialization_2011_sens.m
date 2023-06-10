% Initialize the data

function [] = mf_initialization_2011() % Default value is coeff=1
    filename = 'y_Sobie_ss_2011.mat';
    
    if constants.ind == 0
        Po = 0;
    else
        Po = 5/28;
    end
    
    c_ca_lumen = 1000;
    c_ca_ss = 0.1;
    
    btot = [24 47 900];  
    kp = 1e-3 * [100 115 115];
    km = 1e-3 * [38 100 1000]; 
    b = btot.*(km./kp)./(km./kp+c_ca_ss);
        
    y_Sobie_ss = [Po,c_ca_lumen,c_ca_ss,b];
    % y_Sobie_ss = [Po,c_ca_lumen,c_ca_ss];
    
    
    save(filename,'y_Sobie_ss');
end