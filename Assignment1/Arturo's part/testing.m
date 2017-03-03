%% Explicit eulers
close all; clear all;
[tnList,ynList] = ExplicitEulers(@func,2,25,1);

hold on
    plot(ynList)
    plot(exp(1*tnList)*1)  %e^(lambda*x)*x0
    legend('aprox','real')
hold off

% line width 3
% font size 16

%% Impicit Euler
close all; clear all;
[tnList,ynList] = ImplicitEulers(@func,3.5,100,1);
[tnList1,ynList1] = ExplicitEulers(@func,3.5,100,1);

hold on
    plot(tnList,ynList)
    plot(tnList1,ynList1)
    plot(tnList,exp(1*tnList)*1)  %e^(lambda*x)*x0
    legend('implicit','explicit','real')
hold off