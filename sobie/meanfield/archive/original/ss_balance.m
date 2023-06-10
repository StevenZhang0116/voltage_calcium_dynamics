function [dt_ca_ss] = ss_balance(dhpr_status, j_release, j_efflux, j_buf)
% Equation 1
% Input: Different fluxes determined in different channels / diffusions
%        dhpr_status: whether the dhpr channel is open
% Output: Dynamics of Ca2+ concentration in the subspace
    
    i_dhpr = -0.5 * 10^-12; % DHPR single-channel current A
    v_ss = 1.0 * 10^-13 * 10^-6; % Subspace volume L
    F = 96480; % Faraday constant C * mol^-1
    
    if dhpr_status == 1 % Equation 2
        j_dhpr = -1 * i_dhpr / (2 * F * v_ss); % Fluxes through DHPR channel
    else
        j_dhpr = 0; % if it is closed, the flux is set as 0
    end
        
    dt_ca_ss = j_release + j_dhpr + j_efflux + j_buf;
    
end