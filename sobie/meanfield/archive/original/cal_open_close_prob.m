function [k_open, k_close] = cal_open_close_prob(cal, oxy, n_tot)
% Equation 13 & 14
% Calculate the open/close probablity of the channel based on current
% concentration
% Input: Calcium channel and the total amount of channels

    k_open = oxy + 1 / n_tot;
    k_close = cal + 1 / n_tot;
end

