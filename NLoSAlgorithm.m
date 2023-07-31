% Author: Paulo Francisco
% Objective: NLoS estimator
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

function [error, pos_m, pos_s, u, v, k]=NLoSAlgorithm(b,toas,aods,aoas,m,c)
    
    [u,v]=determinePts(aods,aoas,toas*c,b);
    pos_m=lineIntersect3D(u,v);
    u=u';
    v=v';
    pos_m=pos_m';
    [pos_s,k]=estimateSCs (b,pos_m,aoas,aods,toas);
    
    error=norm(pos_m-m);
    
end

function [y, k]=estimateSCs(EB,EM,aoas,aods,toas)
    ns=length(toas);
    pi=[EB';EM'];
    for i=1:ns
        s(:,i)=determineU(aods(i,1),aods(i,2),toas(i),EB);
        k(:,i)=determineK(aoas(i,1),aoas(i,2),toas(i),EM);
        pf=[s(:,i)';k(:,i)'];
        y(:,i)=lineIntersect3D(pi,pf);
    end    
end

function [s,r]=determinePts(aods,aoas,toas,EB)
    for i=1:size(aods,1)
        s(i,:)=determineU(aods(i,1),aods(i,2),toas(i),EB);
        r(i,:)=determineV(aoas(i,1),aoas(i,2),toas(i),EB); 
    end
end

function p=determineU(aa,ae,toa,EB)
    m=[sin(ae)*cos(aa);
        sin(ae)*sin(aa);
        cos(ae)];
    p=toa*m+EB;
end

function p=determineV(aa,ae,toa,EB)
    m=[sin(pi-ae)*cos(aa-pi);
        sin(pi-ae)*sin(aa-pi);
        cos(pi-ae)];
    p=toa*m+EB;
end

function p=determineK(aa,ae,toa,EM)
    m=[sin(pi-ae)*cos(aa-pi);
        sin(pi-ae)*sin(aa-pi);
        cos(pi-ae)];
    p=EM-toa*m;
end