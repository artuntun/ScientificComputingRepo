function [value,isterminal,direction] = myEventsFcn(t,y)
% when value is equal to zero, an event is triggered.
% set isterminal to 1 to stop the solver at the first event, or 0 to
% get all the events.
%  direction=0 if all zeros are to be computed (the default), +1 if
%  only zeros where the event function is increasing, and -1 if only
%  zeros where the event function is decreasing.
if y(1)<10e-5 || y(2)<10e-5
    value = 0;
else
    value = 1;
end
isterminal = 1;        % Stop the integration
direction = 0;         % All direction

end

