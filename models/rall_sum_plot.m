%% implementation of rall model

L=2;
del=1/16;
%ND

[~,x]=cable_rall(L,del,3,1);
v1=x(:,1);
[~,x]=cable_rall(L,del,3,2);
v2=x(:,1);
[~,x]=cable_rall(L,del,3,3);
v3=x(:,1);
[t,x]=cable_rall(L,del,3,4);
v4=x(:,1);
v_sum=v1+v2+v3+v4;

figure(1)
plot(t,v1,'LineWidth',2,'Color', uint8([50 50 50]))
hold on
plot(t,v2,'LineWidth',2,'Color', uint8([100 100 100]))
plot(t,v3,'LineWidth',2,'Color', uint8([150 150 150]))
plot(t,v4,'LineWidth',2,'Color', uint8([200 200 200]))
%plot(t,v_sum,'k--','LineWidth',2)
xlabel('T/\tau')
ylabel('V')
set(gcf,'unit','normalized','position',[0,0.1,0.35,0.4])
xlim([0,5])
set(gca,'Fontsize',25)
title(['L=',num2str(L),', delay=',num2str(del)])
[~,x]=cable_rall(L,del,0,0);
v_ND=x(:,1);
plot(t,v_ND,'b','LineWidth',2)
[~,x]=cable_rall(L,del,1,0);
v_PD=x(:,1);
plot(t,v_PD,'r','LineWidth',2)
%title('Non-prefered direction (proximal to distal)')
%legend('2&3','4&5','6&7','7&8','sum','ND','PD')
legend('2&3','4&5','6&7','7&8','ND','PD')



%PD
figure(2)
[~,x]=cable_rall(L,del,3,1);
v1=x(:,1);
[~,x]=cable_rall(L,del,3,2);
v2=x(:,1);
[~,x]=cable_rall(L,del,3,3);
v3=x(:,1);
[t,x]=cable_rall(L,del,3,4);
v4=x(:,1);
plot(t,v1,'LineWidth',2)
hold on
plot(t,v2,'LineWidth',2)
plot(t,v3,'LineWidth',2)
plot(t,v4,'LineWidth',2)
plot(t,v_sum,'k--','LineWidth',2)
xlabel('T/\tau')
ylabel('V')
set(gcf,'unit','normalized','position',[0,0.1,0.4,0.4])
set(gca,'Fontsize',30)
title(['L=',num2str(L),', delay=',num2str(del)])
[t,x]=cable_rall(L,del,1,0);
plot(t,x,'k--','LineWidth',2)
v_PD=x(:,1);
%title('Prefered direction (distal to proximal)')
legend('2&3','4&5','6&7','7&8','sum','PD sequence')


figure(3)
plot(t,v1,'LineWidth',2)
hold on
plot(t,v2,'LineWidth',2)
plot(t,v3,'LineWidth',2)
plot(t,v4,'LineWidth',2)
plot(t,v_sum,'LineWidth',2)
plot(t,v_ND,'LineWidth',2)
ylim([0,0.2])
xlabel('T/\tau')
ylabel('V')
set(gcf,'unit','normalized','position',[0,0.1,0.5,0.5])
set(gca,'Fontsize',30)
title('Non-prefered direction (proximal to distal)')
legend('2&3','4&5','6&7','7&8','sum','ND sequence')

% figure(4)
% plot(t,v1,'LineWidth',2)
% hold on
% plot(t,v2,'LineWidth',2)
% plot(t,v3,'LineWidth',2)
% plot(t,v4,'LineWidth',2)
% plot(t,v_sum,'LineWidth',2)
% plot(t,v_ND,'LineWidth',2)
% ylim([0,0.2])
% xlabel('T/\tau')
% ylabel('V')
% set(gcf,'unit','normalized','position',[0,0.1,0.5,0.5])
% set(gca,'Fontsize',30)
% title('Non-prefered direction (proximal to distal)')
% legend('2&3','4&5','6&7','7&8','sum','ND sequence')
% %title(['L=',num2str(L),', delay=',num2str(del)])