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
[tnList,ynList] = ImplicitEulers(@func,@Jacob,100,10000,1);
[tnList1,ynList1] = ExplicitEulers(@func,100,10000,1);

hold on
    plot(tnList,ynList)
    plot(tnList1,ynList1)
    plot(tnList,exp(1*tnList)*1)  %e^(lambda*x)*x0
    legend('implicit','explicit','real')
hold off