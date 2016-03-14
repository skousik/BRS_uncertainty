function sim_conf_ex2
    x = linspace(-0.999,0.999,101);
    theta = linspace(-0.5,0.5,101);
    [X,Y,Z]=meshgrid(x,x,theta);
    XYZ = [X(:) Y(:) Z(:)]';
    Checks = zeros(size(XYZ,2),1);
    options = odeset('Events',@stopp);
    parfor i = 1:size(XYZ,2)
        z0 = XYZ(:,i);
            [~,q] = ode45(@dynamics, [0,1], z0, options);
            final_x = q(end,1:2);
            if norm(final_x) <= 0.5
                Checks(i,1) = 1;
            end
%             keyboard
    end
%     keyboard
    Checks = reshape(Checks,size(X));
%     keyboard
    Check_sum = sum(Checks,3);
%     f = (Check_sum == length(theta));
    
    surf(x,x,Check_sum/numel(theta))
    view(2)
    keyboard
%     z0 = [1e-1,1e-1];
%     [Step.t,Step.q] = ode45(@dynamics, [0,50], z0);
%     plot(Step.q(:,1),Step.q(:,2))
end

function dz = dynamics(~,z)
    xA = z(1:2);
    theta = z(3);
    dx = [ -2*xA(2); 0.8*xA(1)+(9+5*theta)*(xA(1)^2-0.21)*xA(2)];
    dtheta = 0;
    dz = [dx; dtheta];
end

function [value,isterminal,direction] = stopp(~,z)
x = z(1:2);
theta = z(3);
isterminal = 1;
direction = 0;
value = +(max(abs(x))>=1);
end