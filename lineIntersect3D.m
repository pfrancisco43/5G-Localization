% Author: Paulo Francisco
% Objective: Find intersection point of lines in 3D space, in the least squares sense.
% Syntax:
%       [P_intersect] = lineIntersect3D(PA,PB)
% Inputs:
%       PA:     Nx3-matrix containing starting point of N lines
%       APB:    Nx3-matrix containing end point of N lines
%
% Outputs:
%       P_intersect: Best intersection point of the N lines, in least squares sense.

%Algorithm (1) from the paper
function [P_intersect] = lineIntersect3D(PA,PB)    
    Si = PB - PA; %N lines described as vectors
    L=size(PA,1);
    for i=1:L
        ss=sqrt(sum(Si(i,:).^2));
        for j=1:3
            Si(i,j)=Si(i,j)/ss;
        end
    end
    for i=1:3     
        ss=0;
        for j=1:3
            if i==j
                H(i,j)=sum(Si(:,j).^2-1);
                ss=ss+PA(:,j).*(Si(:,j).^2-1);
            else
                H(i,j)=sum(Si(:,i).*Si(:,j));
                ss=ss+PA(:,j).*(Si(:,i).*Si(:,j));
            end            
        end      
        c(i,1)=sum(ss);
    end    
    P_intersect = (H\c)';
end