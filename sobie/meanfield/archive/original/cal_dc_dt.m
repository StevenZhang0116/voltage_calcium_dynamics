function [dc_dt] = cal_dc_dt(cal, k_close, k_open)
% Additional equation
% Nonlinear ODE to calculate change of calcium concentration
% Input: Probablity of channel open/close and current concentration of
%        calcium
% Output: Change of calcium concentration

    dc_dt = k_close * (1 - cal) - k_open * cal;
end

