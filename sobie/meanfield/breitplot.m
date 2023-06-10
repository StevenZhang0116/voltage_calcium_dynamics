% plot the experimental data of spine mentioned in the second paper of
% Breit's group

close all

files = dir('breitdata/*.xls');
aa = cell(1,length(files));
name = cell(1,length(files));

for i = 1:length(files)
    file = files(i);
    aa{i} = readmatrix(file.name);
    name{i} = file.name; 
end

figure()
for i=1:length(aa)
    hold on
    spinedata = aa{i};
    spinedata = sortrows(spinedata);
    xx = spinedata(:,1);
    yy = spinedata(:,2);

    plotyy = smooth(smooth(yy));

%     % iterate through diff interpolation degree
%     degreeset = 1:1:100;
%     errorset = [];
%     for jj=1:length(degreeset)
%         degg = degreeset(jj);
%         xx = xx(~isnan(xx));
%         yy = yy(~isnan(yy));
%         [pp,SS] = polyfit(xx,yy,degg);
%         [~,err] = polyval(pp,xx,SS);
%         errorset(end+1) = sqrt(sum(err.^2));
%     end
% 
%     [~,indd] = min(errorset);
%     optdeg = degreeset(indd);
%     optdeg
%     pp = polyfit(xx,yy,optdeg);
%     plotyy = polyval(pp,xx);

    plot(xx,plotyy,'LineWidth',2);
    xlim([0,20])
    ylim([0 20])
end

bestdata = load('bestdata.mat');
plot(bestdata.t,bestdata.cal_ss,'LineStyle','-.','LineWidth',3);

% sobierosado = load('sobie-rosado-result.mat');
% plot(sobierosado.t,sobierosado.ssv,'LineStyle','-.','LineWidth',3)

legend({'spine2','spine3','spine5','spine6','spine 7',...
    'simulation from Sobie 2022','sobie-rosado'}, 'FontSize',12)
xlabel('Time (ms)','FontSize',16)
ylabel('Concentration [\mu M]','FontSize',16)
title('Comparison','FontSize',16)









