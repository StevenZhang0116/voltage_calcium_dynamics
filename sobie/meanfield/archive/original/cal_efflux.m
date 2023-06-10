function [j_efflux] = cal_efflux(ca_ss)
% Equation 6
% Input: ca_ss: Defined in sobie.m
% Output: Efflux of Ca2+ from the subspace via diffusion to the myoplasm
    
    ca_myo = 0.1 * 10^-6;  % Bulk myoplasmic Ca2+ concentration M 
    tao_efflux = 7.0 * 10^-7; % SS efllux time constant s
    j_efflux = (ca_myo - ca_ss) / tao_efflux;
    
end
