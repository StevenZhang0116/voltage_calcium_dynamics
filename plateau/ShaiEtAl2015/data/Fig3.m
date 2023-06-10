close all
clear all

vs_c = load('O_vs_C.dat');
vd_c = load('O_vd_C.dat');

vs_noca = load('O_vs_noca.dat');
vd_noca = load('O_vd_noca.dat');

vs_c_onlyS = load('O_vs_C_onlyS.dat');
vd_c_onlyS = load('O_vd_C_onlyS.dat');

vs_c_onlyD = load('O_vs_C_onlyD.dat');
vd_c_onlyD = load('O_vd_C_onlyD.dat');

vs_c_DC = load('O_vs_C_DC');
vd_c_DC = load('O_vd_C_DC');

vs_noca_DC = load('O_vs_noca_DC.dat');
vd_noca_DC = load('O_vd_noca_DC.dat');

time = 0:0.05:1000;

time = time(1:length(vs_c));

subplot(2,4,1)
plot(time,vs_c,'k')
hold on
plot(time,vd_c,'r')
xlim([0 200]);ylim([-100 50])
title('Control')

subplot(2,4,2)
plot(time,vs_noca,'k')
hold on
plot(time,vd_noca,'r')
xlim([0 200]);ylim([-100 50])
title('GABA')

subplot(2,4,3)
plot(time,vs_c_onlyS,'k')
hold on
plot(time,vd_c_onlyS,'r')
xlim([0 200]);ylim([-100 50])
title('Only Soma')

subplot(2,4,4)
plot(time,vs_c_onlyD,'k')
hold on
plot(time,vd_c_onlyD,'r')
xlim([0 200]);ylim([-100 50])
title('Only Dend')

subplot(2,4,5)
plot(time,vs_c_DC,'k')
hold on
plot(time,vd_c_DC,'r')
xlim([0 200]);ylim([-100 50])
title('Hyperpolarize')

subplot(2,4,6)
plot(time,vs_noca_DC,'k')
hold on
plot(time,vd_noca_DC,'r')
xlim([0 200]);ylim([-100 50])
title('Hyperpolarize and GABA')