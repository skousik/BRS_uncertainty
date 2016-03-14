% Implementing Dubin's car model for paper
clear all
clc
degree = 6;
T = 1;
freeFinalTime = 0;

nominal_theta = pi/2;
thetabnds = [nominal_theta - 0.01;nominal_theta + 0.01];
v = 1.75; % speed in m/s

% define Cartesian target set
% Xtar = [-2;2];
% Ytar = [-2;2];
% Thetatar = [-pi/32; pi/32];

% change to z-coordinates
% z1tar = Thetatar
% z2tar = [Xtar(1)*cos(Thetatar(1)) + Ytar(1)*sin(Thetatar(2));...
%     Xtar(2)*cos(Thetatar(2)) + Ytar(2)*sin(Thetatar(2))]
% z3tar = [Xtar(2)*sin(Thetatar(1)) - Ytar(2)*cos(Thetatar(2));...
%     Xtar(2)*sin(Thetatar(2)) - Ytar(1)*cos(Thetatar(2))]

% z1tar = Thetatar;
% z2tar = [-1;1];
% z3tar = [-1;1];


% abnds = [v*cos(thetabnds);v]
% bbnds = [v*sin(-thetabnds);v*sin(thetabnds)]

% z1bnds = [-pi/32; pi/32]
z1bnds = [-2;2]
z2bnds = [-2;2]

% define optimization variables
t = msspoly('t',1);

z{1} = msspoly('za',2);
theta{1} = msspoly('pa',1);

zA = z{1};
pA = theta{1};
f{1} = T*[v - zA(2)*pA; zA(1)*pA; 0];

hX{1,1} = [-(zA(1) - z1bnds(1))*(zA(1) - z1bnds(2)); -(zA(2) - z2bnds(1))*(zA(2)-z2bnds(2))];

hTheta{1,1} = [-(theta{1} - thetabnds(1))*(theta{1} - thetabnds(2))];

hXT{1,1} = [0.25 - zA'*zA];

dl{1} = boxMoments(z{1}, [z1bnds(1);z2bnds(1)], [z1bnds(2);z2bnds(2)]);
dtheta{1} = boxMoments(pA, thetabnds(1), thetabnds(2));

% optimize
[w,sol] = solve_BRS_conf(t,z,theta,f,hX,hXT,hTheta,dtheta,dl,degree);
subs(w{1},[z{1};pA],[-1;0;nominal_theta])
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
z1 = linspace(z1bnds(1),z1bnds(2),51);
z2 = linspace(z2bnds(1),z2bnds(2),51);
% z4 = linspace(-thetabnds,thetabnds,11);

[X,Y] = meshgrid(z1,z2);
XY = [X(:) Y(:)]';
P = [1 0; 0 -1];
hold on
Wrestr = msubs(w{1},[z{1};theta{1}],[XY; nominal_theta*ones(size(X(:)))']);
Wrestr_CARTESIAN = msubs(w{1},[z{1};theta{1}],[P*XY; nominal_theta*ones(size(X(:)))']);
Wrestr = reshape(Wrestr,size(X));
Wrestr_CARTESIAN = reshape(Wrestr_CARTESIAN,size(X));
% contour(z1,z2,Wrestr,[1,1],'b')
contour(z1,z2,Wrestr_CARTESIAN,[1,1],'b');
TargetSet = msubs(hXT{1,1},z{1},XY);
TargetSet = reshape(TargetSet,size(X));
contour(z1,z2,TargetSet,[0,0],'r');

%% Plot trajectories
TSPAN = [0 T];
Y0 = [-1; 0];
omega = nominal_theta;

N = 6; %generate <= N^2 points within ROA

% Get point within ROA

[z2_ROA, z1_ROA] = find(Wrestr>=1);
min_z1 = z1(min(z1_ROA(3:end-11)));
max_z1 = z1(max(z1_ROA(11:end-2)));
min_z2 = z2(min(z2_ROA(11:end-11)));
max_z2 = z2(max(z2_ROA(11:end-11)));

z1plot = linspace(min_z1,max_z1,N);
z2plot = linspace(min_z2,max_z2,N);
hold on

% Z1 = repmat(z1plot,N,1);
% Z2 = repmat(z2plot,N,1);
% new_Wrestr = msubs(w{1},[z{1};theta{1}], [Z1(:) Z2(:) nominal_theta*ones(size(Z1(:)))]');
% z1vals = Z1(new_Wrestr >= 1);
% z2vals = Z2(new_Wrestr >= 1);

for i = 1:length(z1plot)
    for j = 1:length(z2plot)
    Y0 = [z1plot(i); z2plot(j)];
    [TOUT, YOUT] = ode45(@(t,y) [v - y(2)*omega; y(1)*omega], TSPAN, Y0);
    
    heading = TOUT*nominal_theta;
    for k = 1:length(TOUT)
        Q = [cos(heading(k)) sin(heading(k)); sin(heading(k)) -cos(heading(k))];
        YOUT(k,:) = (Q\YOUT(k,:)')';
    end
    
    plot(YOUT(1,1),YOUT(1,2),'bo')
    plot(YOUT(:,1),YOUT(:,2),'Color',[0.5;0.5;0.5])
    plot(YOUT(end,1),YOUT(end,2),'mx')
    end
end