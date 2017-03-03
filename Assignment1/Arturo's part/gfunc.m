function F = gfunc(t,yn1,yn,h)
%Return F(y) = yn1 - yn - h*f(yn1,tn1)
f = feval(func,t,y);
F =yn1-yn-h*f;
end
