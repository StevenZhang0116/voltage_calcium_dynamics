clear;

tot_result = tot_equation2;

I_DHPR_func = @(t)0.0-0.5*(t>0.02)*(t<0.0205);
ryr_sobie_test = @(t,y) tot_result.test_ss_model(y,I_DHPR_func(t));

load('y_Sobie_ss.mat')
y0 = y_Sobie_ss;

opts = odeset('RelTol',1e-9,'AbsTol',1e10,'MaxStep',1e-4);
[t,y] = ode45(ryr_sobie_test,[0,0.1],y0,opts);

j_ryr_lst = 0*t;
for i=1:length(t)
    [dy,j_ryr] = tot_result.test_ss_model(y(i,:),I_DHPR_func(t(i)));
    j_ryr_lst(i) = j_ryr;
end

figure(1);
plot(t*1e3,y(:,1));