% Author: Paulo Francisco
% Objective: Verify if exists LoS path
% Syntax:
%       [haveLos]=verifyLos(AoDs, AoAs, Tol) 
% Inputs:
%       AoDs - all angles of departure
%       AoAs - all angles of arrival
%       Tol - tolerance value for considering the existence of a LoS path
%
% Outputs:
%       haveLos - 0 (indicates that only NLoS paths exist)
%                 1 (indicates that a LoS path exists)

function [haveLos]=verifyLos(AoDs, AoAs, Tol)    
    aod=AoDs(1,:);
    aoa=AoAs(1,:);
    
    if aoa(1)<=0.5*pi
        opAz=aoa(1)+pi; 
    else
        opAz=aoa(1)-pi; 
    end
    opEl=pi-aoa(2); 

    haveLos = abs(aod(1)-opAz) <= Tol && abs(aod(2)-opEl) <= Tol; 
end
