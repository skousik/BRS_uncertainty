function endpoints = findDubinsTrajectoryEndpoints(Nw, T, thetaDist)
    endpoints = zeros(2,Nw) ;
    idx = 1 ;

    TSPAN = [0 T]; % time span
    Y0 = [0;0] ; % initial condition
    D = thetaDist ;
    for omega = 1:length(D)
        [TOUT, YOUT] = ode45(@(t,y) [v - y(2)*D(omega); y(1)*D(omega)], TSPAN, Y0);

        heading = TOUT*D(omega);
        for k = 1:length(TOUT)
            Q = [cos(heading(k)) sin(heading(k)); sin(heading(k)) -cos(heading(k))];
            YOUT(k,:) = (Q\YOUT(k,:)')';
        end

        endpoints(1,idx) = YOUT(end,1) ;
        endpoints(2,idx) = YOUT(end,2) ;

        idx = idx + 1 ;
    end
end