function [tnList,ynList] = ImplicitEulers1(func,tspan,N,Y0, Niter)
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
sizeY = size(Y0);
Ndim = sizeY(1);
sizeTspan = size(tspan);
Ninit = sizeTspan(2);  % We can be given only the end time, then the begining is 0

ynList = zeros(Ndim,N+1);
tnList = zeros(Ndim,N+1);

%% Initialization
if (Ninit == 1)
    tbegin = 0;
    tend = tspan;
elseif(Ninit == 2) % [tbegin tend]
    tbegin = tspan(1,1);
    tend = tspan(1,2);
end
dt = (tend - tbegin)/N;

tnList(1) = tbegin;
ynList(:,1) = Y0;

%% Loop

for k = 1:N
    yn1Est = ynList(:,k);
    for it = 1:Niter
%         if (it == 1) || (it == Niter)
%             yn1Est
%         end
        f = feval(func,tnList(k) + dt,yn1Est);
        yn1Est = ynList(:,k)+dt*f;
    end
    ynList(:,k+1) = yn1Est;
    tnList(k+1) = tnList(k)+dt;
end