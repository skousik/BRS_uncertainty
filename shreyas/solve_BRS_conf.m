function [w,sol] = solve_BRS_conf(t,x,theta,f,hX,hXT,hTheta,mu_theta,dl,degree)
% number of hybrid modes
nModes = length( x );

% define the time constraint
T = 1;
hT  = t*(T-t);

% define the program
prog = spotsosprog;
prog = prog.withIndeterminate( t );

for i = 1:nModes
    % make x{i} indeterminate
    prog = prog.withIndeterminate( x{ i } );
    prog = prog.withIndeterminate( theta{ i } );
    
    % create the v{i} variable
    vmonom{ i } = monomials( [ t; x{ i }; theta{i}], 0:degree );
    [ prog, v{ i }, ~ ] = prog.newFreePoly( vmonom{ i } );     
  
    % create the w{i} variable
    wmonom{ i } = monomials( [x{ i };theta{i}], 0:degree );
    [ prog, w{ i }, wcoeff{ i } ] = prog.newFreePoly( wmonom{ i } );

    % creating the variables that will be used later
    v0{ i } = subs( v{ i }, t, 0 );
    vT{ i } = subs( v{ i }, t, T );
    dvdt{ i } = diff( v{ i }, t );
    dvdx{ i } = diff( v{ i }, [x{ i }; theta{i}] );
    Lfv{ i } = dvdt{ i } + dvdx{ i } * f{ i };     
end

% creating the constraints and cost function
obj = 0;
for i = 1:nModes
        % for perfect
        prog = sosOnK( prog, -Lfv{ i }, [ t; x{ i };theta{i}], [ hT; hX{ i };hTheta{i}], degree );
        [prog,m01] = sosOnK( prog, w{ i } - v0{i} -1, x{ i }, [hX{ i }; hTheta{i}], degree );
        prog = sosOnK( prog, w{ i }, x{ i }, hX{ i }, degree );
        
        % mu_T1
        [prog,mT] = sosOnK( prog, vT{i}, [x{ i };theta{i}], [hXT{ i };hTheta{i}], degree );

        obj = obj + mu_theta{i}( dl{ i }(w{i}) );
end

% set options
options = spot_sdp_default_options();
options.verbose = 1;

% actually solve the problem
sol = prog.minimize( obj, @spot_mosek_sos, options );
for i = 1:nModes
    w{ i } = sol.eval( w{ i } );
end
% keyboard