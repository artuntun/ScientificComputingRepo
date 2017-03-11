 function [tnList,ynList] = ImplicitEulersAdaptiveStep(func,Jacob,tspan,N,Y0,abstol,reltol)
%
% This function solves a general first-order Initial Value Problem
% of the form
%                Ydot = F(y,t),  y(tbegin) = Y0
%
% using Implicit Euler Method in N steps (constant step size)
%
% INPUT:
%    func   : a function handle to function F(Y,t)
%    Jacob  : a function handle to function dF(Y,t)/dY
%    tspan  : a 1x2 array of the form [tbegin tend]
%    N      : used to calculate the intial time step size
%    Y0     : initialvalue(s)s
%    abstol : maxmium absolute error perimtted
%    reltol : maxmium relative error perimtted
%

sizeY = size(Y0);
Ndim = sizeY(1);
sizeTspan = size(tspan);
Ninit = sizeTspan(2);

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
ynList(:,1)=Y0;

tol = 10e-5; % Newtons method tolerance

%%error estimation and control paramters
epstol = 0.8;
facmin = 0.1;
facmax = 5;

%% Loop
k=1;
while tnList(k) < tend
    
    %In order to compute until tend
    if (tnList(k)+h>tend)
        h = tend-tnList(k);
    end
    
    f = feval(func,tnList(k),ynList(:,k));
    
    acceptedStep = 0;
    while ~acceptedStep
        
        %y(n+1) in one step
        y_guess = ynList(:,k)+h*f; %Forwards euler
        y = newtonODE(func,Jacob,y_guess,tol,tnList(k)+h,ynList(:,k),h);
        
        %y(n+1) in two step
        hm = h/2;
        ym_guess = ynList(:,k)+hm*f;
        ym = newtonODE(func,Jacob,ym_guess,tol,tnList(k)+hm,ynList(:,k),hm);
        f = feval(func,tnList(k)+hm,ym);
        y_hat_guess = ynList(:,k)+hm*f;
        y_hat = newtonODE(func,Jacob,y_hat_guess,tol,tnList(k)+h,ym,hm);
        
        e = abs(y_hat - y);
        r = max(e./max(abstol,abs(y_hat).*reltol));
        if r<=1.0
            acceptedStep = 1;
            ynList(:,k+1) = y_hat;
            tnList(k+1) = tnList(k)+h;
        end
        %hnew = h*sqrt(wanted_error/current_error)
        h = max(facmin, min(sqrt(epstol/r),facmax))*h;
    end
    k = k+1;
end