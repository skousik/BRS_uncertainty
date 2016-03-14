function dzdt = dubins_dzdt( t,z,v,omega )
% creates a function to be used in ODE45
% assuming theta is not a state. so the system becomes a two state system.
dzdt = [v - z(2)*omega; z(1)*omega];

end

