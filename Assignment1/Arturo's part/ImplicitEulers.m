function [tnList,ynList] = ImplicitEulers(func,Jacob,tspan,N,y0)
%
% This function solves a general first-order Initial Value Problem
% of the form
%                u? = f(u,t),  u(tstart) = eta
%
% using Euler?s Method in n steps (constant step size).s
%
% INPUT:
%    func  : a function handle to function f(u,t)
%    tspan : a 1x2 array of the form [tstart tend]
%    N     : total number of steps in tspan
%    y0   : initialvalue(s)s
%    param : parameters to be passed to func
%
ynList = zeros(1,N+1);
tnList = zeros(1,N+1);
dt = tspan/N;
tol = 10e-4;
ynList(1)=y0;

for k = 1:N
    f = feval(func,tnList(k),ynList(:,k));
    y_guess = ynList(:,k)+dt*f; %Forwards euler
    yn1 = newtonODE(y_guess,tol,tnList(k)+dt,ynList(k),dt,func,Jacob);
    f = feval(func,tnList(k)+dt,yn1);
    ynList(k+1) = ynList(k)+dt*f;
    tnList(k+1) = tnList(k)+dt;
end