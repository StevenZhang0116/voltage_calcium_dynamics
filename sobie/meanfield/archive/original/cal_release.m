function [j_release] = cal_release(num_open, ca_lumen, ca_ss) 
% Equation 3-4
% Input: ca_lumen: Ca2+ concentration in the local, JSR lumen
%        ca_ss: Ca2+ concentration in a restricted ubsarcolemmal space
%        num_open: Num of RyR channels that are open
%        num_open + num_close = num_total, which is 20/50/100 normally
% Output: j_ryr: Flux through a single, open RyR channel

    d_ryr = 4000; % Ca2+ diffusion through open RyR parameter 1/s
    j_ryr = d_ryr * (ca_lumen - ca_ss); 
    
    j_release = num_open * j_ryr; % change the num_open by open probablity
                                  % described by ODE
    
end
