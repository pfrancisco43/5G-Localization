% Author: Paulo Francisco
% Objective: Show Scenario of LoS Localization
% Syntax:
%       showLoSScenario(b,m,pos)
% Inputs:
%       b, m, pos - BS, true MS and estimated MS coordinates

function showLoSScenario(b,m,pos)
    figure
    plot3(b(1),b(2),b(3),'b^');
    hold on
    plot3(m(1),m(2),m(3),'rx');
    plot3(pos(1),pos(2),pos(3),'k.')

    line([b(1) b(1)],[b(2) b(2)],[b(3) 0],'LineStyle','--','Color','b')
    line([b(1) pos(1)],[b(2) pos(2)],[b(3) pos(3)],'Color','b')
    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')

    title('Scenario for LoS Estimation')
    legend('BS','MS','Estimated MS','Location','north')
    view(45,25)
end