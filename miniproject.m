%% Challange 21.1
clear all; close all;
k = 4;
ro = 0.8;
y0 = [0.005 1-0.005 0];
options = odeset('Events',@myEventsFcn);

[t,y] = ode23(@(t,y) infect(y,ro,k),[0 800],y0,options);

plot(y)
legend('I','S','R')