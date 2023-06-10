jsr1 = load("2002_jsr.mat").CaJSR_all;
jsr2 = load("2011_jsr.mat").CaJSR_avg;
jsr3 = load("2011_jsr.mat").CaJSRtot_all * 1000;
load("2002_plottime.mat")
% 
% plottime = plottime(1:1000);
% jsr1 = jsr1(1:1000);
% jsr2 = jsr2(1:1000);
% jsr3 = jsr3(1:1000);
figure()
hold on
plot(plottime, jsr1)
plot(plottime, jsr2)
plot(plottime, jsr3)
plot(plottime, (1 - exp(-plottime/90)) * 1000,'r','LineWidth',2.5)
lgd = legend('2002','2011','2011\_total','2011\_total asymptote');
lgd.FontSize = 14;
% title("Number of Open Channels")
xlabel("Time(ms)")
ylabel("JSR [Ca^{2+}] (\muM)")
title("JSR [Ca^{2+}]")
hold off

open1 = load("2002_numopen.mat").Nopen_all_2002;
open2 = load("2011_numopen.mat").Nopen_all_2011;
% load("2002_plottime.mat")
% % CaJSR_all = CaJSR_all/1000;
% 
% figure()
% hold on
% plot(plottime, open1)
% plot(plottime, open2)
% lgd = legend('2002','2011');
% lgd.FontSize = 14;
% title("Number of Open Channels")
% xlabel("Time(ms)")
% ylabel("Number of Open RyR Channels (RyR Open Probablity)")
% hold off

a = load("2002_releaseflux.mat");
b = load("2011_releaseflux.mat");
load("2002_plottime.mat")
% CaJSR_all = CaJSR_all/1000;

figure()
hold on
plot(plottime, a.Irel_all)
plot(plottime, b.Irel_all)
lgd = legend('2002','2011');
lgd.FontSize = 14;
title("SR Ca^{2+} Release Flux")
xlabel("Time(ms)")
ylabel("SR Ca^{2+} Release Flux")
hold off

%% Calculation
cnt1 = [50,0];
cnt2 = [50,0];
for i = 1:1000
    if i > 50 && a.Irel_all(i) ~= 0
        cnt1(1) = i;
        break
    end
end
for i = 1:1000
    if i > 50 && b.Irel_all(i) ~= 0
        cnt2(1) = i;
        break
    end
end

for i = 1:1000
    if i > 55 && a.Irel_all(i) == 0
        cnt1(2) = i;
        break
    end
end
for i = 1:1000
    if i > 55 && b.Irel_all(i) == 0
        cnt2(2) = i;
        break
    end
end
    
%% integral of release flux
result1 = 0; 
result2 = 0;
for i = cnt1(1):cnt1(2)
    if i == cnt1(1) || i == cnt1(2)
        result1 = result1 + a.Irel_all(i);
    else
        result1 = result1 + 2 * a.Irel_all(i);
    end
end

for i = cnt2(1):cnt2(2)
    if i == cnt2(1) || i == cnt2(2)
        result2 = result2 + b.Irel_all(i);
    else
        result2 = result2 + 2 * b.Irel_all(i);
    end
end

disp(result1/result2)

%% integral of jsr
temp1=0; temp2=0;
for i=1:length(jsr1)
    if(jsr1(i)==min(jsr1))
        temp1 = i;
        break
    end
end

for i=1:length(jsr2)
    if(jsr2(i)==min(jsr2))
        temp2 = i;
        break
    end
end

ijsr1 = 1000-jsr1;
ijsr2 = 1000-jsr2;

% temp1=cnt1(2);
% temp2=cnt2(2);
int1=0; int2=0;
for i=52:temp1
    if i==cnt1(1) || i==temp1
        int1 = int1 + ijsr1(i);
    else
        int1 = int1 + 2*ijsr1(i);
    end
end

for i=52:temp2
    if i==cnt2(1) || i==temp2
        int2 = int2 + ijsr2(i);
    else
        int2 = int2 + 2*ijsr2(i);
    end
end

disp(int1/int2)


%% Calculate current through a single RyR
cnt = 0;
for i = 1:1000
    if jsr2(i) ~= 1000
        cnt = i;
        break;
    end
end

totalc = b.Irel_all(cnt);
totalc/open2(cnt)


