% Author: Paulo Francisco
% Objective: get antenna array (UCA)
% Syntax:
%       array = getArray(N)
% Inputs:
%       N:  Number of antennas
%
% Outputs:
%       array - antenna array
%       
function array = getArray(N)
    radius = 0.4/sind(180/N);
    rx = radius*cosd(360*(0:N-1).'/N);
    ry = radius*sind(360*(0:N-1).'/N);
    r = [rx, ry];
    array = r;
end
