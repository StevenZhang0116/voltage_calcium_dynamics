%% Draw the graph for raw/relative sensitivity analysis (for single parameter)

function [timelst] = raw_sens_graph(name)
    data = load(name).all_result;
    xdata = data(:,1);
    % Figure 1 
    figure()
    set(gca, 'YScale', 'log')
    hold on
    timelst = [];
    for i=2:4
        plot(xdata,data(:,i),'LineWidth',2);
        [m,indx] = max(data(:,i));
        time = xdata(indx);
        timelst(end+1) = time;
    end
    legend('Open Prob','Ca2+ Lumen', 'Release flux')
    title(name)
    xlabel('k')
    hold off
%     figname1 = "fig1"+name+".png";
%     saveas(gcf,figname1);
    % Figure 2
    figure()
    set(gca, 'YScale', 'log')
    hold on
    for i=5:7
        plot(xdata,data(:,i),'LineWidth',2);
        [m,indx] = max(data(:,i));
        time = xdata(indx);
        timelst(end+1) = time;
    end
    legend('Open Prob','Ca2+ Lumen', 'Release flux')
    title(name)
    xlabel('k')
    hold off
%     figname2 = "fig2"+name+".png";
%     saveas(gcf,figname2);
    % Figure 3
%     norm_coeff = abs(log(xdata));
%     avg1 = (data(:,2)+data(:,5))/2;
%     avg2 = (data(:,3)+data(:,6))/2;
%     avg3 = (data(:,4)+data(:,7))/2;
%     figure(3)
%     set(gca, 'YScale', 'log')
%     hold on
%     plot(xdata,avg1./norm_coeff,'LineWidth',2);
%     plot(xdata,avg2./norm_coeff,'LineWidth',2);
%     plot(xdata,avg3./norm_coeff,'LineWidth',2);
%     legend('turb ratio - Open Prob','turb ratio - Ca2+ Lumen', 'turb ratio - release flux')
%     title('turb ratio '+name);
%     hold off
%     disp(norm_coeff)
end