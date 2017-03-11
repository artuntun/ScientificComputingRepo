function [tnList,ynList] = RK4AdaptiveStep(func,tspan,N,Y0,abstol,reltol)
%
% This function solves a general first-order Initial Value Problem
% of the form
%                u? = f(u,t),  u(tstart) = eta
%
% using Euler?s Method in n steps (constant step size).
%
% INPUT:
%    func  : a function handle to function f(u,t)
%    tspan : a 1x2 array of the form [tstart tend]
%    N     : total number of steps in tspan
%    y0   : initialvalue(s)
%    param : parameters to be passed to func
%
sizeTspan = size(tspan);
Ninit = sizeTspan(2);  % We can be given only the end time, then the begining is 0

%% Initialization
if (Ninit == 1)
    tbegin = 0;
    tend = tspan;
elseif(Ninit == 2) % [tbegin tend]
    tbegin = tspan(1,1);
    tend = tspan(1,2);
end
h = (tend - tbegin)/N;

ynList = [];
tnList = [];
tnList(1) = tbegin;
ynList(:,1) = Y0;

%%error estimation and control paramters
epstol = 0.8;
facmin = 0.1;
facmax = 5;
kpow = 0.2;

%% Loop
k=1;
while tnList(k) < tend
    %In order to compute until tend
    if (tnList(k)+h>tend)
        h = tend-tnList(k);
    end
    
    acceptedStep = 0;
    while ~acceptedStep
        
        %y(n+1) in one step
        y=RungeKuttaStep(func,tnList(k),ynList(:,k),h);
        
        %y(n+1) in two step
        hm = h/2;
        ym=RungeKuttaStep(func,tnList(k),ynList(:,k),hm);
        y_hat=RungeKuttaStep(func,tnList(k)+hm,ym,hm);
        
        e = abs(y_hat - y);
        r = max(e./max(abstol,abs(y_hat).*reltol));
        if r<=1.0
            acceptedStep = 1;
            ynList(:,k+1) = y_hat;
            tnList(k+1) = tnList(k)+h;
        end
        %hnew = h*sqrt(wanted_error/current_error)
        h = max(facmin, min(sqrt(epstol/r)^kpow,facmax))*h;
    end
    k = k+1;
end