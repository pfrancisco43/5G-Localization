% Author: Paulo Francisco
% Objective: get the response vector
% Syntax:
%       resp=getResponse(array,angs)
% Inputs:
%       array: Antenna array
%       angs: pair of the azimuth and elevation angles
%
% Outputs:
%       resp - response vector
%      
function resp=getResponse(array,angs)    
    ang=angs;
    N=length(array);
    rho = [cos(ang(1))*cos(ang(2)); sin(ang(1))*cos(ang(2))];  
    resp = exp(-1j*pi*array*rho)/sqrt(N);
end