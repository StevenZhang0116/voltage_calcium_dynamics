%% Relative sensitivity analysis
% check1 = standard, check2 = variation

function result = rela_analysis(check1,check2)    
    indexlst = [2,3,8]; 
    
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
        result(end+1) = val;
    end

end


