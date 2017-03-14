function [F] = IVP2(t,Y,params)
%IVP number 2 
    A = [0 1;-1 0];  
    F = A*Y;
end

