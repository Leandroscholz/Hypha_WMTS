function [value,isterminal,direction] = myEvent(t,y)
global Deltax

% Locate the time when nth reactor length is equal to twice length of
% previous tanks 
% and stop integration.
value = y(end)-2*Deltax;     % detect length = 2*Deltax
isterminal = 1;   % stop the integration
direction = 0;    % negative direction