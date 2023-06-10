function [bi_dt] = cal_bi_dt(ca_ss)
% Additional Equation
% Not sure if I understand this equation correctly...
% If so, it resembles to equation (5) very much

    b_cam_tot = 24.0 * 10^-6; % Total calmodulin concentration M
    k_cam_on = 100 * 10^-6; % Calmodulin Ca2+ on-rate constant M^-1 * s^-1
    k_cam_off = 38.0; % Calmodulin Ca2+ off-rate constant
    
    b_sr_tot = 47.0 * 10^-6; % Total SR membrane concentration M
    k_sr_on = 115 * 10^-6; % SR membrane Ca2+ on-rate constant M^-1 * s^-1
    k_sr_off = 100.0; % SR membrane Ca2+ off-rate constant
    
    b_sl_tot = 1124.0 * 10^-6; % Total calmodulin concentration M
    k_sl_on = 115 * 10^-6; % Calmodulin Ca2+ on-rate constant M^-1 * s^-1
    k_sl_off = 1000.0; % Calmodulin Ca2+ off-rate constant
    
    % These three items are never numerically defined in the paper
    % Meaning as unbounded buffer concentration
    % First order ODE related to ca_ss over time
    b_cam_unbounded = 0.8 * b_cam_tot;
    b_sr_unbounded = 0.8 * b_sr_tot;
    b_sl_unbounded = 0.8 * b_sl_tot;
    
    bi_dt = - k_cam_on * ca_ss * b_cam_unbounded + k_cam_off ...
        * (b_cam_tot - b_cam_unbounded) - k_sr_on * ca_ss ...
        * b_sr_unbounded + k_sr_off * (b_sr_tot - b_sr_unbounded)...
        - k_sl_on * ca_ss * b_sl_unbounded + k_sl_off ...
        * (b_sl_tot - b_sl_unbounded);
end

