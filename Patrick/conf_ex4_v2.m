% Implementing Dubin's car model for paper
clear all
clc
degree = 8;
T = 1;
freeFinalTime = 0;

thetabnds = .001;
v = 1; % speed in m/s

% define Cartesian target set
Xtar = [-1;1];
Ytar = [-1;1];
Thetatar = [-pi/32; pi/32];

% change to z-coordinates
% z1tar = Thetatar
% z2tar = [Xtar(1)*cos(Thetatar(1)) + Ytar(1)*sin(Thetatar(2));...
%     Xtar(2)*cos(Thetatar(2)) + Ytar(2)*sin(Thetatar(2))]
% z3tar = [Xtar(2)*sin(Thetatar(1)) - Ytar(2)*cos(Thetatar(2));...
%     Xtar(2)*sin(Thetatar(2)) - Ytar(1)*cos(Thetatar(2))]

z1tar = Thetatar;
z2tar = [-1;1];
z3tar = [-1;1];


% abnds = [v*cos(thetabnds);v]
% bbnds = [v*sin(-thetabnds);v*sin(thetabnds)]

z1bnds = [-pi/32; pi/32]
z2bnds = [-3;3]
z3bnds = [-3;3]

% define optimization variables
t = msspoly('t',1);

z{1} = msspoly('za',3);
theta{1} = msspoly('pa',1);

zA = z{1};
pA = theta{1};
f{1} = T*[pA; v - zA(3)*pA; zA(2)*pA; 0];

hX{1,1} = [-(zA(1) - z1bnds(1))*(zA(1) - z1bnds(2)); -(zA(2) - z2bnds(1))*(zA(2)-z2bnds(2));...
    -(zA(3) - z3bnds(1))*(zA(3) - z3bnds(2))];

hTheta{1,1} = [(theta{1}+thetabnds)*(thetabnds-theta{1})];

hXT{1,1} = [1 - zA(2:end)'*zA(2:end)];

dl{1} = boxMoments(z{1}, [z1bnds(1);z2bnds(1);z3bnds(1)], [z1bnds(2);z2bnds(2);z3bnds(2)]);
dtheta{1} = boxMoments(pA, -thetabnds, thetabnds);

% optimize
[w,sol] = solve_BRS_conf(t,z,theta,f,hX,hXT,hTheta,dtheta,dl,degree);
subs(w{1},[z{1};pA],[0;0;0;0])
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

%% Discretizing space

% z1 = linspace(z1bnds(1),z1bnds(2),21);
z2 = linspace(z2bnds(1),z2bnds(2),51);
z3 = linspace(z3bnds(1),z3bnds(2),51);
% z4 = linspace(-thetabnds,thetabnds,11);

[X,Y] = meshgrid(z2,z3);
XY = [X(:) Y(:)]';
subplot(211)
hold on
Wrestr = msubs(w{1},[z{1};theta{1}],[zeros(size(X(:)))'; XY; zeros(size(X(:)))']);
Wrestr = reshape(Wrestr,size(X));
contour(z2,z3,Wrestr,[1,1]);
TargetSet = msubs(hXT{1,1},[z{1};theta{1}],[zeros(size(X(:)))'; XY; zeros(size(X(:)))']);
TargetSet = reshape(TargetSet,size(X));
contour(z2,z3,TargetSet,[0,0],'r');
