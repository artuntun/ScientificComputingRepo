function [ gradY ] = VanDerPol( t, Y, args )
%VANDELPOL Summary of this function goes here
%   Detailed explanation goes here
    mu = 100;
    % Syntax: xdot = PreyPredator(t,x,a,b)
    gradY = zeros(2,1);
    gradY(1) = Y(2);
    gradY(2) = mu*(1-power(Y(1),2))*Y(2) - Y(1);

end

