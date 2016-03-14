function uncertaintyDist = randomSampleUncertaintyDist(theta, thetaBounds, ...
                                                     shape, N_theta)
    if iscell(theta)
        theta = theta{1} ;
    end
    
    thBoundLt = thetaBounds(1) ;
    thBoundRt = thetaBounds(2) ;

    % start with uniform distribution
    theta_distribution = linspace(thBoundLt, thBoundRt, N_theta) ;
    
    % values of shape over a range of thetas
    ftheta_vec = msubs(shape,theta,theta_distribution)';

    uncertaintyDist = zeros(N_theta,1);
    % 
    randos = rand([N_theta,1]);
    deltaTheta = thBoundRt - thBoundLt ;
    probs = cumsum(ftheta_vec*(deltaTheta/(N_theta-1)));

    for i = 1:length(randos)
        uncertaintyDist(i) = theta_distribution(find(probs>=randos(i),1)) + ...
                             (deltaTheta/(N_theta-1))/2;
    end
end