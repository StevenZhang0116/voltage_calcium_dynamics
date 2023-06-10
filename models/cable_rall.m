%% reproducing rall's result of compartment model 

function [t,x]=cable_rall(L,del,direction,loc)

    N=10; % number of total compartments
    tspan=0:0.1:5;
    x0=zeros(N,1);  % initial condition
    [t,x]=ode23(@(t,x) MOL_ODE(t,x),tspan,x0);
    % plot the membrane potential time courses
    % figure(1)
    % plot(t,x(:,1),'LineWidth',2)
    % hold on
    % xlabel('T/\tau')
    % ylabel('V')
    % set(gcf,'unit','normalized','position',[0,0.1,0.4,0.4])
    % set(gca,'Fontsize',30)
    % title(['L=',num2str(L),', delay=',num2str(del)])
    
    % figure
    % hold on
    % for i=(N+2):(2*N+2)
    %     plot(t,x(:,i))
    % end

    function dx=MOL_ODE(t,x)
        N=length(x);
        %L=5;
        delta_x=L/N;
        %del=0.25;
        if direction==0
            ton=[0;0;del;del;2*del;2*del;3*del;3*del];
        elseif direction==1
            ton=[3*del;3*del;2*del;2*del;del;del;0;0];
        elseif direction==2
            ton=[100;100;100;100;100;100;100;100];
            ton(loc*2)=(loc-1)*del;
            ton(loc*2-1)=(loc-1)*del;
            %ton=[100;0;0;100;100;100;100;100;100;100];
            %ton=[100;100;100;del;del;100;100;100;100;100];
            %ton=[100;100;100;100;100;2*del;2*del;100;100;100];
            %ton=[100;100;100;100;100;100;100;3*del;3*del;100];
        elseif direction==3
            ton=[100;100;100;100;100;100;100;100];
            ton(loc*2)=(4-loc)*del;
            ton(loc*2-1)=(4-loc)*del;
            %ton=[100;3*del;3*del;100;100;100;100;100;100;100];
            %ton=[100;100;100;2*del;2*del;100;100;100;100;100];
            %ton=[100;100;100;100;100;del;del;100;100;100];
            %ton=[100;100;100;100;100;100;100;0;0;100];
        end
        toff=ton+0.25;
        E=heaviside(t-ton).*heaviside(toff-t);
        dx=x*0;
        dx(1)=-x(1)+(2*x(2)-2*x(1))/(delta_x^2)-0*(x(1)-1);
        dx(2:N-1)=-x(2:N-1)+(x(1:N-2)-2*x(2:N-1)+x(3:N))/(delta_x^2)-E.*(x(2:N-1)-1);
        dx(N)=-x(N)+(2*x(N-1)-2*x(N))/(delta_x^2)-0*(x(N)-1);
    end

end
