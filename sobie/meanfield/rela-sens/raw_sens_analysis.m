%% Brute force sensitivity analysis
%% 2011 Model using 2002 Parameter

% 1. Ratio of Integral
% 2. Relative distance when t is the same


function [result,sens] = raw_sens_analysis(name,coeff)
    sens = 0;

    check1 = load('meanfield_result_2011_standard.mat').result;
    check2 = load(name).result;
    
    indexlst = [2,3,8]; %% respectively of open prob, lumen concentration, release flux
    
    time1 = check1(:,1);
    time2 = check2(:,1);

    result = [];
    
    for i = 1:length(indexlst)
        compare1 = check1(:,indexlst(i));
        compare2 = check2(:,indexlst(i));
    
        %% integral diff
        integral1 = 0; integral2 = 0;
        % trapezoidal rule to calculate riemann sum
        if(indexlst(i) == 3) % Rotate the graph
            compare1 = 1000 - compare1;
            compare2 = 1000 - compare2;
        end
    
        for j = 2:length(time1)
            integral1 = integral1 + (time1(j)-time1(j-1))*(compare1(j)+compare1(j-1))/2;
        end
        for j = 2:length(time2)
            integral2 = integral2 + (time2(j)-time2(j-1))*(compare2(j)+compare2(j-1))/2;
        end
%         figure()
%         hold on
%         plot(time1,compare1,'LineWidth',1)
%         plot(time2,compare2,'LineWidth',1);
%         legend('original','modified')
%         hold off
    
        result(end+1) = abs(integral1-integral2)/integral1;
    end

    sum = 0; % check whether stable
    
    for i = 1:length(indexlst)
        compare1 = check1(:,indexlst(i));
        compare2 = check2(:,indexlst(i));
    
        %% point by point diff
        val1 = max(compare1);
        val2 = max(compare2);
        err_lst=[];
    
%         if(indexlst(i) == 3) 
%             compare1 = 1000 - compare1;
%             compare2 = 1000 - compare2;
%         end
    
        threshold = 0.001; %% ignore the values that are very small (very unstable)
    
        for j=1:min(length(time1),length(time2))
            if(compare1(j)>threshold*val1 && compare2(j)>threshold*val2 && compare1(j) > 0 && compare2(j) > 0)
                err = abs((compare1(j)-compare2(j)))/compare1(j); 
                err_lst(end+1) = err;
            end
        end
        val = mean(err_lst);
        if (val < 5*10e-3) sum = sum + 1; end
        result(end+1) = val;
    end
%     disp(result)
    if (sum == 3)
        sens = coeff;
    end
end


