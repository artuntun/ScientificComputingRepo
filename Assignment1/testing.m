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
[tnList,ynList] = ImplicitEulers(@func,2,25,1);

hold on
    plot(ynList)
    plot(exp(1*tnList)*1)  %e^(lambda*x)*x0
    legend('aprox','real')
hold off