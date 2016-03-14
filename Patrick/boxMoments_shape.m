function m = boxMoments_shape(x,f,a,b)
%
% m = boxMoments_shape(x,a,b)
% f is the RN derivative of the mu wrt lambda
% Produces m(p), p being an msspoly, which integrates out the variables
%   x over the box a(i) <= x(i) <= b(i).
    function l = moments(p)
        l = p*f;
        for i = 1:length(x)
            l = subs(integral(l,x(i)),x(i),b(i)) ...
                - subs(integral(l,x(i)),x(i),a(i));
        end
    end
    m = @moments;
end