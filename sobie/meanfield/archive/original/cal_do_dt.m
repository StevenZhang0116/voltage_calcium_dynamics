function [do_dt] = cal_do_dt(cal, k_close, k_open)
% Additional equation
% Nonlinear ODE to calculate change of oxygen concentration
% Input: Probablity of channel open/close and current concentration of
%        calcium
% Output: Change of oxygen concentration

    do_dt = k_open * cal - k_close * (1 - cal);
end

