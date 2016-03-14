function sim_conf_ex1
    x = linspace(0,1,101);
    theta = linspace(-.5,0.5,101);
    [X,Y]=meshgrid(x,theta);
    XY = [X(:) Y(:)]';
    Checks = zeros(size(XY,2),1);
    for i = 1:size(XY,2)
        z0 = XY(:,i);
            [Step(1).t,Step(1).q] = ode45(@dynamics, [0,1], z0);
            final_x = Step(1).q(end,1);
            if final_x <= 0.6 && final_x >=0
                Checks(i,1) = 1;
            end
    end
%     keyboard
    Checks = reshape(Checks,size(X));
    Check_sum = sum(Checks,1);
    f = (Check_sum == length(theta));
    
    figure
    subplot(311)
    plot(x,f,'*',[-1,-1],[0,1],'k:',[1,1],[0,1],'k:')
    ylim([0,1.5])
    xlim([-1.1,1.1])
    subplot(312)
    surf(x,theta,reshape(Checks,size(X)))
    view(2)
    subplot(313)
    plot(x,Check_sum/numel(theta))
    keyboard
end

function dz = dynamics(~,z)
    x = z(1);
    theta = z(2);
    dx = -x^2*theta-0.5*x^3;
    dtheta = 0;
    dz = [dx; dtheta];
end