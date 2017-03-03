function dydt = Jacob(t,y,h,func)
%Retrun the function dy/dt=f(t,y(t))
f = feval(func,t,y);
dydt=1-h*f;
end

