function [ xn1 ] = newtonRootFinding(x0,tol)
%UNTITLED20 Summary of this function goes here
%   Detailed explanation goes here
xn = x0;
maxIt = 100;
c = 2;
b = 1;
[f,fprime]=f1(x0);
err = 1000;
it = 0;

while err>=tol && it < maxIt
    it = it +1;
    [f,fprime]=f1(xn);
    xn1 = xn - (f/fprime);
    [f,fprime]=f1(xn1);
    err = f;
    xn = xn1;
end

end

