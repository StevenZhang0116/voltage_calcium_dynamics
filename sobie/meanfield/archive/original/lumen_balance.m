function [dt_ca_lumen] = lumen_balance(j_ryr, ca_lumen)
% Equation 7-9
% Describe the balance equation for JSR lumenal Ca2+ concentration 
% Input: j_ryr: Defined in sobie.m 
%        ca_lumen: Defined in sobie.m
% Output: Dynamics of JSR lumenal Ca2+ concentration
    
    csq_total = 10.0 * 10^-3; % Total calsequestrin concentration M 
    k_csq  = 0.8 * 10^-3; % Calsequestrin Ca2+ dissociation constant M 
    csq_unbounded = caq_total; % This term is never defined explicitly in
                               % in Sobie's paper, but used in the
                               % equation, so set by default for now
    beta_jsr = 1 / (1 + (csq_total * k_csq) / (k_csq + csq_unbounded)^2);
    % Buffering in the JSR by calsequestrin
    
    v_ss = 1.0 * 10^-13 * 10^-6; % Subspace volume L
    v_jsr = 1.0 * 10^-11 * 10^-6; % JSR volume L 
    
    ca_nsr = 1.0 * 10^-3; % NSR Ca2+ concentration M 
    tau_refill = 0.01; % JSR refilling time constant s
    j_refill = (ca_nsr - ca_lumen) / tau_refill;
    
    dt_ca_lumen = beta_jsr * (-j_ryr * v_ss / v_jsr + j_refill);
    
end
