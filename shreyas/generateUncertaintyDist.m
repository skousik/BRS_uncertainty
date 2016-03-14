function [shape, shape_int] = generateUncertaintyDist(theta, thetaBounds, weights)%, root_locations)

% TO DO:
% - check total degree of poly
% - generate roots either according to user-input locations or autospaced
% - generate shape function according to those roots

    if length(thetaBounds) == 2
        thetaBounds = sort(thetaBounds) ;
        thetaBoundLt = thetaBounds(1) ;
        thetaBoundRt = thetaBounds(2) ;
    else
        error('thetabounds must be a 2-element vector of doubles')
    end    
    
    if iscell(theta)
        theta = theta{1} ;
    end
    
%     if nargin < 4
%         % if root_locations isn't provided, evenly space the roots, based on
%         % the number of weights, between the theta bounds
%         root_locations = linspace(thetaBoundLt,thetaBoundRt,length(weights)) ;
%     end

    l_wt = weights(1) ;
    r_wt = weights(end) ;

    if length(weights) == 3
        c_wt = weights(2) ;
    else
        c_wt = 0 ;
    end

        
    % shape determines distribution of uncertainty
    shape = -(theta)^(c_wt)*...
             (theta-thetaBoundLt)^(l_wt)*...
             (theta-thetaBoundRt)^(r_wt);

     % function to integrate shape
    shape_int = boxMoments_shape(theta,shape,thetaBoundLt,thetaBoundRt);
    
    volume = double(shape_int(msspoly(1))) ;
    
    % create final outputs
    shape = shape/volume ;
    shape_int = boxMoments_shape(theta,shape,thetaBoundLt,thetaBoundRt);
end