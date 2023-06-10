%% Analayze the parameter combinations output

result = load('result.mat');
set_para = result.stable_lst;
indicator = result.all_result;
all_para = load("paras.mat").all_parameters;

five_percentage = length(set_para) / length(all_para); 

one_para = zeros(1,11);
one_ind = zeros(1,6);

threshold = 0.05;

% the order of parameters: 
% alpha, dryr, nryr, csq, efflux, refill, vjsr, vss, bsl, kcsq, kmax
for i = 2: length(set_para)
    this_para = set_para(i,:);
    this_ind = indicator(i,:);
    if(this_ind(4) < threshold && this_ind(5) < threshold && this_ind(6) < threshold)
        one_ind(end+1,:) = this_ind;
        one_para(end+1,:) = this_para;
    end
end

length(one_para)-1

% one_para
