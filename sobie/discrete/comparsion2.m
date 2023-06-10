ryr_lst = load('N_RyR_lst.mat').N_RyR_lst;
endtime1 = load('2002_total_endtime.mat').total_endtime;
fullopentime1 = load('2002_total_fullopentime.mat').total_fullopentime;
endtime2 = load('2011_total_endtime.mat').total_endtime;
fullopentime2 = load('2011_total_fullopentime.mat').total_fullopentime;

figure()
hold on
plot(ryr_lst,endtime1)
plot(ryr_lst,fullopentime1)
plot(ryr_lst,endtime2)
plot(ryr_lst,fullopentime2)
xlabel("Total RyR Channels")
ylabel("Time(ms)")
legend("2002 Total End Time","2002 Full Open Time","2011 Total End Time",...
    "2011 Full Open Time")
lgd.FontSize = 14;
% title("RyR Channel Total End Time & Full Open Time for 2002 & 2011 Models")
hold off