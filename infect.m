function dydt = infect(y,ro,k)
%Retrun the infected population model
%   Detailed explanation goes here

dydt=[ro*y(1)*y(2)-y(1)/k; -ro*y(1)*y(2); y(1)/k];
end

