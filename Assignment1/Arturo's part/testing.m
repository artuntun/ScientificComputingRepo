% line width 3
% font size 16
%% Impicit Euler
close all; clear all;
N = 100;
tspan = [0 10];
[tnList,ynList] = ImplicitEulers(@func,@Jacob,tspan,N,1);
[tnList1,ynList1] = ExplicitEulers(@func,tspan,N,1);
[tnList2,ynList2] = ImplicitTrapezoid(@func,@Jacob,tspan,N,1);

hold on
    plot(tnList,ynList,'LineWidth',2)
    plot(tnList1,ynList1)
    plot(tnList2,ynList2,'--','LineWidth',2)
    plot(tnList,exp(-1*tnList)*1)  %e^(lambda*x)*x0
    legend('implicit','explicit','real')
hold off
%% Implict Euler for multidimensional for Predator-Prey
close all; clear all;
tspan = [0 200];
Y0=[0.5; 0.5];
N = 10000;
[tnList,ynList] = ImplicitEulers(@DepPrey,@JacobDepPrey,tspan,N,Y0);
[tnList1,ynList1] = ExplicitEulers(@DepPrey,tspan,N,Y0);

hold on
    plot(ynList(1,:),ynList(2,:),'LineWidth',1)
    plot(ynList1(1,:),ynList1(2,:),'LineWidth',1)
    leg = legend('Implict','Explicit');
    set(leg,'FontSize',16);
    grid
hold off
%% Explicit Eulers with adaptive step
close all; clear all;
abstol = 1e-3;
reltol = 1e-4;
tspan = [0 50];
Y0 = 1;
% [tnList2,ynList2] = ImplicitEulers(@func,@Jacob,20,100,1);
[tnList1,ynList1] = ExplicitEulersAdaptiveStep(@func,tspan,1000,...
    Y0,abstol,reltol);
[tnList,ynList] = ExplicitEulers(@func,tspan,1000,Y0);

hold on
%     plot(tnList2,ynList2)
    plot(tnList1,ynList1,'LineWidth',2)
    plot(tnList,ynList,'LineWidth',2)
    plot(tnList,exp(1*tnList)*1,'LineWidth',2)  %e^(lambda*x)*x0
    leg = legend('adaptive','fixstep','real');
    set(leg,'FontSize',16);
    grid
hold off
%% Implicit Eulers: fix step vs adaptive step
close all; clear all;
abstol = 1e-3;
reltol = 1e-4;
tspan = [0 50];
N = 1000;
Y0 = 1;
[tnList1,ynList1] = ImplicitEulersAdaptiveStep(@func,@Jacob,tspan,...
    N,Y0,abstol,reltol);
[tnList,ynList] = ImplicitEulers(@func,@Jacob,tspan,N*2,Y0);

hold on
    plot(tnList1,ynList1,'LineWidth',2)
    plot(tnList,ynList,'LineWidth',2)
    plot(tnList,exp(1*tnList)*1,'LineWidth',2)  %e^(lambda*x)*x0
    leg = legend('Adaptive','Fix-step','real');
    set(leg,'FontSize',16);
    grid
hold off

%% Testin Classic Runge Kutta

close all; clear all;
abstol = 1e-3;
reltol = 1e-4;
tspan = [0 50];
N = 1000;
Y0 = 1;
[tnList,ynList] = RK4AdaptiveStep(@func,tspan,N,Y0,abstol,reltol);

hold on
    plot(tnList,ynList,'LineWidth',2)
    plot(tnList,exp(1*tnList)*1,'LineWidth',2)  %e^(lambda*x)*x0
    leg = legend('RK','real');
    set(leg,'FontSize',16);
    grid
hold off