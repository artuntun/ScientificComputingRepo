function [tnList,ynList] = ImplicitEuler(func,tspan,N,y0)
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
%    y0   : initialvalue(s)
%    param : parameters to be passed to func
%

ynList = zeros(1,N+1);
tnList = zeros(1,N+1);
dt = tspan/N;

ynList(1)=y0;

for k = 1:N
    yn1 = newtonODE(ynList(k)+0.05,tol,tnList(k)+dt,ynList(k),dt);
    ynList(k+1) = yn1;
    tnList(k+1) = tnList(k)+dt;
end