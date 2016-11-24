function [I] = Current(t)
%CURRENT Summary of this function goes here
%   Detailed explanation goes here

%Change the current after 1/2 of the time period
if (t < 12.5)
I = 0.1;
% The current value
else
I = 0.3; 


end

