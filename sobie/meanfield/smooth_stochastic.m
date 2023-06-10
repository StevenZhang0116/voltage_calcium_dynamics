a = load('2011_numopen.mat').Nopen_all_2011;
b = load('2011_jsr.mat').CaJSR_avg;
c = load('2011_releaseflux.mat').Irel_all;
time = load('2002_plottime.mat').plottime;

mf_result_1 = load('meanfield_result_2011_1.mat').result;
mf_result_0 = load('meanfield_result_2011_0.mat').result;

figure()
hold on
plot(time, a/constants.n_ryr);
plot(mf_result_1(:,1),mf_result_1(:,2));
plot(mf_result_0(:,1),mf_result_0(:,2));
xlim([0,100])
xlabel('Time(ms)');
ylabel('Open Probablity');
legend('200 Stochastic Independent Trails average',...
    'Mean-field Representation DHPR initialization',...
    'Mean-field Representation Original initialization');
hold off

figure()
hold on
plot(time,b);
plot(mf_result_1(:,1),mf_result_1(:,3));
plot(mf_result_0(:,1),mf_result_0(:,3));
xlim([0,100])
xlabel('Time(ms)');
ylabel('JSR Ca^{2+} concentration')
legend('200 Stochastic Independent Trails average',...
    'Mean-field Representation DHPR initialization',...
    'Mean-field Representation Original initialization');
hold off

figure()
hold on
plot(time,c);
plot(mf_result_1(:,1),mf_result_1(:,8));
plot(mf_result_0(:,1),mf_result_0(:,8));
xlim([0,100])
xlabel('Time(ms)');
ylabel('SR Ca^{2+} Release Flux')
legend('200 Stochastic Independent Trails average',...
    'Mean-field Representation DHPR initialization',...
    'Mean-field Representation Original initialization');
hold off

max(c)
max(mf_result_1(:,8))
max(mf_result_0(:,8))