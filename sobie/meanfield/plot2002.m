%--------------------------------------------------------------------------
% Plot the comparison result for smoothing DHPR current
%
% Steven Zhang, Courant Institute
% Updated July 2022
%--------------------------------------------------------------------------

close all

path = './res/';
S = dir(fullfile(path,'**','*.mat'));
names = {S.name};

% param
int = 0.0001;
xx = 3:int:6;

% generate color
colorlst = [];
for i = 1:length(names)
    colorlst(end+1,:) = [rand,rand,rand] ;
end

timet = [];
figure()
hold on
for i = 1:length(names)
    matname = char(names(i));
    loader = load([path,matname]);

    aa = loader.t > 10;
    for j=1:length(aa) 
        if aa(j) == 1 
            ind1 = j; 
            break; 
        end 
    end
    bb = loader.t > 4;
    for j=1:length(bb) 
        if bb(j) == 1 
            ind2 = j; 
            break 
        end 
    end
    timet(end+1,1) = ind1;
    timet(end,2) = ind2;

    plot(loader.t(ind2:ind1),loader.cal_ss(ind2:ind1),'LineWidth',3, ...
        'Color',colorlst(i,:))
end
hold off
legend(names,'Location','northeast')
xlim([4 10])
title('SR calcium concentration comparsion')
xlabel('time (ms)')
ylabel('concentration (\mu M)')

figure()
hold on
for i = 1:length(names)
    matname = char(names(i));
    loader = load([path,matname]);
    ind1 = timet(i,1);
    ind2 = timet(i,2);
    plot(loader.t(ind2:ind1),loader.prob(ind2:ind1),'LineWidth',3,...
        'Color',colorlst(i,:))
end
hold off
legend(names,'Location','northeast')
xlim([4 10])
title('open prob comparsion')
xlabel('time (ms)')
ylabel('open probablity')

figure()
hold on
xxx = 4:0.001:10;
for i = 1:length(names)
    matname = char(names(i));
    loader = load([path,matname]);
    plot(xxx,loader.I_DHPR_func(xxx),'LineWidth',3,'Color',colorlst(i,:))
end
hold off
legend(names,'Location','northeast')
xlim([4 6.5])
title('DHPR injection')
xlabel('time (ms)')
ylabel('Injected current (nA)')