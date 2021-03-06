close all
clear all
clc

mu = 100;
t0 = 0;
tfinal = max(300);
tspan = [t0 tfinal];

method = 'ESDIRK23';
flag = 0;
isenspar = 1;
x0 = [2; 1];
h0 = 0.001;
absTol = 1e-6;
relTol = 1e-6;
fun  = 'VanderPolFun';
jac  = 'VanderPolJac';

[Tout_ESDIRK23,Xout_ESDIRK23,Gout_ESDIRK23,Eout_ESDIRK23, info,stats] = ESDIRK23(fun,jac,[t0 tfinal],x0,h0,absTol,relTol,method,mu);

nStep = 1:1:info.nStep;


[Tout_RK3,Xout_RK3,Eout_RK3] = RK3(@VanDelPol,tspan,x0,h0, mu);  % @DepPrey @VanDelPol

% method = 'DOPRI54';
% solver = ERKSolverErrorEstimationParameters(method);
% [Tout_RK4,Xout_RK4,Eout_RK4] = ERKSolverErrorEstimation(@VanDelPol,tspan,Y0,h,solver,mu);

lw= 3
figure()
subplot(3,1,1);  
hold on
    plot(Xout_RK3(:,2), Xout_RK3(:,1), 'color',rand(1,3), 'LineWidth',lw)
    plot(Xout_ESDIRK23(:,2), Xout_ESDIRK23(:,1), 'color',rand(1,3), 'LineWidth',lw)
    legend('RK3','ESDIRK23')
    ylabel('x_2','FontSize',12,'FontWeight','bold')
    xlabel('x_1','FontSize',12,'FontWeight','bold')
hold off

subplot(3,1,2); 
hold on
    plot(Tout_RK3, Xout_RK3(:,1), 'color',rand(1,3), 'LineWidth',lw)
    plot(Tout_ESDIRK23, Xout_ESDIRK23(:,1), 'color',rand(1,3), 'LineWidth',lw)
    legend('RK3','ESDIRK23')
    ylabel('x_1','FontSize',12,'FontWeight','bold')
    xlabel('time','FontSize',12,'FontWeight','bold')
hold off

subplot(3,1,3);  
hold on
    plot(Tout_RK3, Xout_RK3(:,2), 'color',rand(1,3), 'LineWidth',lw)
    plot(Tout_ESDIRK23, Xout_ESDIRK23(:,2), 'color',rand(1,3), 'LineWidth',lw)
    legend('RK3','ESDIRK23')
    ylabel('x_2','FontSize',12,'FontWeight','bold')
    xlabel('time','FontSize',12,'FontWeight','bold')
hold off
 
%% Individual figure

figure()
subplot(3,1,1);  
hold on

    plot(Xout_ESDIRK23(:,2), Xout_ESDIRK23(:,1), 'color',rand(1,3), 'LineWidth',lw)
    legend('ESDIRK23')
    ylabel('x_2','FontSize',12,'FontWeight','bold')
    xlabel('x_1','FontSize',12,'FontWeight','bold')
hold off

subplot(3,1,2); 
hold on

    plot(Tout_ESDIRK23, Xout_ESDIRK23(:,1), 'color',rand(1,3), 'LineWidth',lw)
    legend('RK3','ESDIRK23')
    ylabel('x_1','FontSize',12,'FontWeight','bold')
    xlabel('time','FontSize',12,'FontWeight','bold')
hold off

subplot(3,1,3);  
hold on

    plot(Tout_ESDIRK23, Xout_ESDIRK23(:,2), 'color',rand(1,3), 'LineWidth',lw)
    legend('ESDIRK23')
    ylabel('x_2','FontSize',12,'FontWeight','bold')
    xlabel('time','FontSize',12,'FontWeight','bold')
hold off
 

%% Error figure
figure();

subplot(3,1,1); title('Estimated Local Error 1');
hold on
    plot(Tout_RK3, log(abs(Eout_RK3(:,1))), 'color',rand(1,3), 'LineWidth',lw)
    plot(Tout_RK3, log(abs(Eout_ESDIRK23(:,1))), 'color',rand(1,3), 'LineWidth',lw)
    
    legend('RK3','ESDIRK23')
    ylabel('log(|e_1|)','FontSize',12,'FontWeight','bold')
    xlabel('time','FontSize',12,'FontWeight','bold')
    
hold off

subplot(3,1,2); title('Estimated Local Error 2');
hold on
    plot(Tout_RK3, log(abs(Eout_RK3(:,2))), 'color',rand(1,3), 'LineWidth',lw)
    plot(Tout_RK3, log(abs(Eout_RK3(:,2))), 'color',rand(1,3), 'LineWidth',lw)
    legend('RK3','ESDIRK23')
    ylabel('log(|e_2|)','FontSize',12,'FontWeight','bold')
    xlabel('time','FontSize',12,'FontWeight','bold')
    
hold off

subplot(3,1,3); title('Estimated Local Error');
hold on
    plot(Tout_RK3, log(sqrt(Eout_RK3(:,2).^2 + Eout_RK3(:,1).^2 )), 'color',rand(1,3), 'LineWidth',lw)
    plot(Tout_RK3, log(sqrt(Eout_ESDIRK23(:,2).^2 + Eout_ESDIRK23(:,1).^2 )), 'color',rand(1,3), 'LineWidth',lw)

    legend('RK3','ESDIRK23')
    ylabel('log(|Error|)','FontSize',12,'FontWeight','bold')
    xlabel('time','FontSize',12,'FontWeight','bold')
    
hold off


ESDIRKperformance(info, stats);




