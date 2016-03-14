% add the dual folder
addpath('../functions/');
addpath('../');
addpath(genpath('C:\Users\elemsn\Documents\MATLAB\SeDuMi_1_3'))
addpath(genpath('C:\Users\elemsn\Documents\MATLAB\spotless'))
addpath(genpath('C:\Program Files\Mosek\7\toolbox\r2013aom'))
% keyboard
clear all
% close all

degree        = 10;
T             = 1;
freeFinalTime = 0;

t = msspoly( 't', 1 );


x{ 1 }   = msspoly( 'xa', 1 );
theta{1} = msspoly('pa',1);

xA       = x{1};
f{ 1 }   = T *[ -xA^2*theta{1}-0.5*xA^3;0];

hX{1,1} = [(1-xA(1))*(xA(1)+1)];
hXt{1,1} = [(1-xA(1))*(xA(1)+1); (theta{1}+0.5)*(0.5-theta{1})];

hTheta{1,1} = [(theta{1}+0.5)*(0.5-theta{1})];
% hDX{1,1} = [hX{1,1};hTheta{1,1}];
% hDX{1,1} = [xA(1)-1; 0.5-theta{1}; 0.5+theta{1}];
% hDX{1,2} = [-1-xA(1); 0.5-theta{1}; 0.5+theta{1}];
% hDX{1,3} = []


hXT{1,1} = [(xA(1)+0)*(0.6-xA(1))];
% hXTc{1,1} = -hXT{1,1};

dl{ 1 }    = boxMoments( x{ 1 }, -1, 1);
dtheta{ 1 }= boxMoments( theta{1}, -0.5, 0.5);
mu_theta{1} = boxMoments( theta{1}, -.5, 0.5); % you'll have to change this. Consider using the function boxMoments_shape

% dlt{1} = boxMoments([x{1};theta{1}],[-1;-0.5],[1,0.5]);

% sX = cell(1,1);
% R = cell(1,1);

w = solve_BRS_conf(t,x,theta,f,hX,hXT,hTheta,mu_theta,dl,degree);



% [w,v] = solve_BRS_hybrid_inner(t,{[x{1};theta{1}]},f,hXt,hDX,hXTc,sX,R,dlt,degree,freeFinalTime);
%%
x1 = linspace(-1,1,81);
x2 = linspace(-0.5,.5,51);
[X,Y] = meshgrid(x1,x2);
XY = [X(:) Y(:)]';

subplot(211)
hold on
Wrestr = msubs(w{1}^6,[x{1};theta{1}],XY);
contour(x1,x2,reshape(Wrestr,size(X)),[1,1]);
% return
subplot(212)
hold on
Wrestr = 1*msubs(mu_theta{1}(w{1}),x{1},x1);
plot(x1,1-Wrestr/double(mu_theta{1}(msspoly(1))));