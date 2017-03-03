function [ yn1_n1 ] = newtonODE(y0,tol,tn1,yn,dt)
%UNTITLED20 Summary of this function goes here
%   Detailed explanation goes here
yn1_n = y0;
maxIt = 100;
err = 1000;
it = 0;

while err>=tol && it < maxIt
    it = it +1;
    fprime = 1-dt*feval(@func,tn1,yn1_n);
    f = yn1_n-dt*feval(@Jacob,tn1,yn1_n)-yn;
    yn1_n1 = yn1_n - (f/fprime);
    f = yn1_n1-dt*feval(@Jacob,tn1,yn1_n1)-yn;
    err = f;
    yn1_n = yn1_n1;
end

end

