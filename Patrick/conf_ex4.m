% Implementing Dubin's car model for paper
clear all
clc
degree = 8;
T = 1;
freeFinalTime = 0;

thetabnds = .1;
v = 30; % speed in m/s

% define Cartesian target set
Xtar = [29.9;30.1];
Ytar = [-1;1];
Thetatar = [-0.05; 0.05];

% change to z-coordinates
z1tar = Thetatar
z2tar = [Xtar(1)*cos(Thetatar(1)) + Ytar(1)*sin(Thetatar(2));...
    Xtar(2)*cos(Thetatar(2)) + Ytar(2)*sin(Thetatar(2))]
z3tar = [Xtar(2)*sin(Thetatar(1)) - Ytar(2)*cos(Thetatar(2));...
    Xtar(2)*sin(Thetatar(2)) - Ytar(1)*cos(Thetatar(2))]


abnds = [v*cos(thetabnds);v]
bbnds = [v*sin(-thetabnds);v*sin(thetabnds)]

z1bnds = [-thetabnds; thetabnds]
z2bnds = [0;...
    abnds(2)*cos(z1bnds(2)) + bbnds(2)*sin(z1bnds(2))]
z3bnds = [abnds(2)*sin(z1bnds(1)) - bbnds(2)*cos(z1bnds(2));...
    abnds(2)*sin(z1bnds(2)) - bbnds(1)*cos(z1bnds(2))]

% define optimization variables
t = msspoly('t',1);

x{1} = msspoly('xa',3);
theta{1} = msspoly('pa',1);

xA = x{1};
pA = theta{1};
f{1} = T*[pA; v - xA(3)*pA; xA(2)*pA; 0];

hX{1,1} = [-(xA(1) - z1bnds(1))*(xA(1) - z1bnds(2)); -(xA(2))*(xA(2)-z2bnds(2));...
    -(xA(3) - z3bnds(1))*(xA(3) - z3bnds(2))];

hTheta{1,1} = [(theta{1}+thetabnds)*(thetabnds-theta{1})];

hXT{1,1} = [-(xA(1) - z1tar(1))*(xA(1) - z1tar(2)); -(xA(2) - z2tar(1))*(xA(2) - z2tar(2));...
    -(xA(3) - z3tar(1))*(xA(3) - z3tar(2))];

dl{1} = boxMoments(x{1}, [z1bnds(1);z2bnds(1);z3bnds(1)], [z1bnds(2);z2bnds(2);z3bnds(2)]);
dtheta{1} = boxMoments(pA, -thetabnds, thetabnds);

% optimize
[w,sol] = solve_BRS_conf(t,x,theta,f,hX,hXT,hTheta,dtheta,dl,degree);

% %% Time to plot
% 
% % shape = -(pA-0.01)*(pA+0.01)
% shape = 1;
% mu_theta{1} = boxMoments_shape(theta{1},shape,-thetabnds,thetabnds);
% 
% x1 = linspace(-.01,.01,81);
% x2 = linspace(-thetabnds,thetabnds,51);
% [X,Y] = meshgrid(x1,x2);
% XY = [X(:) zeros(size(X(:))) zeros(size(X(:))) Y(:)]';
% 
% subplot(211)
% hold on
% Wrestr = msubs(w{1},[x{1};theta{1}],XY);
% contour(x1,x2,reshape(Wrestr,size(X)),[1,1]);
% % return
% subplot(212)
% hold on
% Wrestr = 1*msubs(mu_theta{1}(w{1}),x{1},[x1,zeros(size(x1)),zeros(size(x1))]);
% plot(x1,Wrestr/double(mu_theta{1}(msspoly(1))));
