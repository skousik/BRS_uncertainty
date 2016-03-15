function [beta, w_plop] = conf_ex4_montecarlo_function(T, v, omegaBounds, ...
                                    R_T, weights, Nx, Ny, Nw, xbds, ybds)
% [beta, w_plop] = conf_ex4_montecarlo_function(T, v, omegaBounds, ...
%                                   R_T, weights, Nx, Ny, Nw, xbds, ybds)
%                               
% Runs a pseudo-Monte-Carlo simulation to generate the BRS of a Dubin's
% car system given yaw rate uncertainty, and returns the alpha-level
% confidence surface as an (Nx x Ny) matrix and the distribution of the
% uncertainty discretized as a (1 x Nw) vector
%
% INPUTS:
%  - T (time scale)
%  - v (vehicle velocity)
%  - omegaBounds (uncertainty bounds as a 2-element vector)
%  - R_T (target set radius)
%  - weights (left and right weights of distribution as a 2-element vec)
%
% OPTIONAL INPUTS:
%  - Nx (number of x points in grid, default 300)
%  - Ny (number of y points in grid, default 300)
%  - Nw (number of omega uncertainty samples, default 300)
%  - xbds (2-element vector of grid x bounds, default [-1, 1])
%  - ybds (2-element vector of grid y bounds, default [-1, 1])
%
% OUTPUTS:
%  - beta (Nx x Ny matrix of alpha-level confidences corresponding to each
%          point on the grid in the xbds/ybds space)
%  - w_plop (1 x Nw vector of the discretized uncertainty shape function)



    % omega distribution
    left_weight = weights(1) ;
    right_weight = weights(end) ;

    if nargin < 6
        % granularity of X, Y, and omega spaces
        Nx = 300 ; % number of x points
        Ny = 300 ; % number of y points
        Nw = 300 ; % number of omegas

        xbds = [-1, 1] ;
        ybds = [-1, 1] ;
    end
    
    x = linspace(xbds(1),xbds(2),Nx) ;
    y = linspace(ybds(1),ybds(2),Ny) ;

    % Create distribution of omega
    % start with uniform distribution
    omegaBdL = omegaBounds(1) ;
    omegaBdR = omegaBounds(end) ;
    omega_distribution = linspace(omegaBdL, omegaBdR, Nw) ;

    % create mss poly to evaluate this thing
    theta{1} = msspoly('pa',1) ; % uncertainty

    l_wt = left_weight ;
    r_wt = right_weight ;
     % shape determines distribution of uncertainty
    shape = -(theta{1}-omegaBdL)^(r_wt)*(theta{1}-omegaBdR)^(l_wt);

     % handle to integrate shape
    mu_theta{1} = boxMoments_shape(theta{1},shape,omegaBdL,omegaBdR);

     % values of shape over a range of thetas
    ftheta_vec = msubs(shape,theta{1},omega_distribution)';

    volume = double(mu_theta{1}(msspoly(1)));
    delta_w = omegaBdR - omegaBdL ;

    num_traj = Nw;
    omegas = zeros(num_traj,1);
    % 
    randos = rand([num_traj,1]);
    probs = cumsum(ftheta_vec*(delta_w/(length(ftheta_vec)-1)))/volume;

    for i = 1:length(randos)
        omegas(i) = omega_distribution(find(probs>=randos(i),1)) + ...
                                       (delta_w/(length(ftheta_vec)-1))/2;
    end

    omega_distribution = omegas ;

    % Determine trajectory endpoints
    endpoints = zeros(2,Nw) ;
    idx = 1 ;

    TSPAN = [0 T]; % time span
    Y0 = [0;0] ; % initial condition
    DIST = omega_distribution ;
    for omega = 1:length(DIST)
        [TOUT, YOUT] = ode45(@(t,y) [v - y(2)*DIST(omega); y(1)*DIST(omega)], TSPAN, Y0);

        heading = TOUT*DIST(omega);
        for k = 1:length(TOUT)
            Q = [cos(heading(k)) sin(heading(k)); sin(heading(k)) -cos(heading(k))];
            YOUT(k,:) = (Q\YOUT(k,:)')';
        end

        endpoints(1,idx) = YOUT(end,1) ;
        endpoints(2,idx) = YOUT(end,2) ;

        idx = idx + 1 ;
    end

    % Monte Carlo sim
    % Create grid of points to fill with tallies via the magical matrix tour
    Z1 = repmat(x,Ny,1,Nw) ;
    Z2 = repmat(y',1,Nx,Nw) ;

    End1 = repmat(reshape(endpoints(1,:),1,1,Nw),Nx,Ny) ;
    End2 = repmat(reshape(endpoints(2,:),1,1,Nw),Nx,Ny) ;

    Z1End = Z1 + End1 ;
    Z2End = Z2 + End2 ;

    mag = sqrt(Z1End.^2 + Z2End.^2) ;

    beta = sum(mag <= R_T, 3) ;
    beta = beta/Nw ;
    
    scale = 1/double(mu_theta{1}(msspoly(1))) ;
    w_plop = ftheta_vec*scale*delta_w ;
end