% Author: Paulo Francisco
% Objective: LoS estimator
% Syntax:
%       [error, pos] = LosAlgorithm(b,toa,aod,m)
% Inputs:
%       b, m - BS and MS coordinates
%       toa - Time of Arrival
%       aod - Angle of Departure
%
% Outputs:
%       error - Eclidian Distance from pos and m
%       pos - MS estimation

function [error,pos]=LoSAlgorithm(b,toa,aod,m,c)    
    x(1)=toa*c*sin(aod(2))*cos(aod(1))+b(1);
    x(2)=toa*c*sin(aod(2))*sin(aod(1))+b(2);
    x(3)=toa*c*cos(aod(2))+b(3);
    pos=x';
    error=norm(pos-m);
end