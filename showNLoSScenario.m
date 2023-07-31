% Author: Paulo Francisco
% Objective: Show Scenario of LoS Localization
% Syntax:
%       showNLoSScenario(b,m,s,pos_m,pos_s,u,v,k)
% Inputs:
%       b, m, s:        BS, MS and SC true cordinates
%       pos_m, pos_s:   MS and SC estimated cordinates
%       u, v, k:        Auxiliar points

function showNLoSScenario(b,m,s,pos_m,pos_s,u,v,k)
    figure
    p1=plot3(b(1),b(2),b(3),'^b','MarkerSize',8);    
    hold on;
    line([b(1) b(1)],[b(2) b(2)],[b(3) -8],'LineStyle','--','Color','b')
    p2=plot3(m(1),m(2),m(3),'rx','MarkerSize',8);  
    p3=plot3(pos_m(1),pos_m(2),pos_m(3),'k.','MarkerSize',8); 
    for i=1:size(s,2)
        p4=plot3(s(1,i),s(2,i),s(3,i),'dr','MarkerSize',8);
        p5=plot3(pos_s(1,i),pos_s(2,i),pos_s(3,i),'ko','MarkerSize',8);   

        %point u
        p6=plot3(u(1,i),u(2,i),u(3,i),'vm');
        p7=line([b(1),u(1,i)],[b(2),u(2,i)],[b(3),u(3,i)],'Color','b');
        %point v
        p8=plot3(v(1,i),v(2,i),v(3,i),'X','Color','#D95319');
        p9=line([v(1,i),u(1,i)],[v(2,i),u(2,i)],[v(3,i),u(3,i)],'Color','b');
        %point k
        %p10=plot3(k(1,i),k(2,i),k(3,i),'X','Color','#D95319');
        %p11=line([k(1,i),pos_m(1)],[k(2,i),pos_m(2)],[k(3,i),pos_m(3)],'Color','red');

        %Propagation path
        %p12=line([b(1),s(1,i)],[b(2),s(2,i)],[b(3),s(3,i)],'LineStyle','--','Color','m','LineWidth',1.5);
        %line([s(1,i),m(1)],[s(2,i),m(2)],[s(3,i),m(3)],'LineStyle','--','Color','m','LineWidth',1.5);
    end
        
    grid on;
    xlabel('x')
    ylabel('y')
    zlabel('z')
     
    title('Scenario for NLoS Estimation')
    leg=legend([p1 p2 p3 p4,p5],{'BS','MS','Estimated MS','SC^k','Estimated SC^k'},'Location','southeast');
    view(75,55);
end