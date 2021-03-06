%% Implementing Dubin's car model for paper
clear all
clc
close all

degree = 10;
T = 1;
freeFinalTime = 0;

nominal_omega = 0; % where the distribution is centered around
omegabnds = [nominal_omega - 3*pi/4;nominal_omega + 3*pi/4]; % bounds of distribution
v = 0.5; % speed in m/s

z1bnds = [-1;1]; % bounds on domain of z1
z2bnds = [-1;1]; % bounds on domain of z2

RT = 0.25;

%% define optimization variables
t = msspoly('t',1); % time

z{1} = msspoly('za',2); % states z1 and z2
theta{1} = msspoly('pa',1); % uncertainty

zA = z{1};
pA = theta{1};
f{1} = T*[v - zA(2)*pA;
              zA(1)*pA;
                     0]; % dynamics

% quadratics define domain of z
hX{1,1} = [-(zA(1) - z1bnds(1))*(zA(1) - z1bnds(2)); -(zA(2) - z2bnds(1))*(zA(2)-z2bnds(2))];

% domain of uncertainty
hTheta{1,1} = [-(theta{1} - omegabnds(1))*(theta{1} - omegabnds(2))];

% target set is a ball around the origin
hXT{1,1} = [RT^2 - zA'*zA];

dl{1} = boxMoments(z{1}, [z1bnds(1);z2bnds(1)], [z1bnds(2);z2bnds(2)]);
dtheta{1} = boxMoments(pA, omegabnds(1), omegabnds(2));

% optimize!!
[w,sol] = solve_BRS_conf(t,z,theta,f,hX,hXT,hTheta,dtheta,dl,degree);

%% Discretizing space

z1 = linspace(z1bnds(1),z1bnds(2),101); % discretize z1
z2 = linspace(z2bnds(1),z2bnds(2),101); % discretize z2
z3 = linspace(omegabnds(1),omegabnds(2),1001); % discretize z3 (domain of omegas)

[X,Y] = meshgrid(z1,z2); % create a grid of z1 and z2
XY = [X(:) Y(:)]'; % this turns the grid into a long vector
P = [1 0; 0 -1]; % transformation matrix from Cartesian to z-coordinates at steering angle = 0 (corresponds to time = 0)

Wrestr_CARTESIAN = msubs(w{1},[z{1};theta{1}],[P*XY; nominal_omega*ones(size(X(:)))']); % this gives you values of w for a grid in X-Y plane
Wrestr_CARTESIAN = reshape(Wrestr_CARTESIAN,size(X)); % reshape vector into grid in X-Y plane

TargetSet = msubs(hXT{1,1},z{1},XY); % define target set over grid
TargetSet = reshape(TargetSet,size(X));

%% Time to integrate dat shit!!
% l_wt = 1;
% r_wt = 1;
% shape = -(theta{1}-omegabnds(1))^r_wt*(theta{1}-omegabnds(2))^l_wt; % shape determines distribution of uncertainty
% mu_theta{1} = boxMoments_shape(theta{1},shape,omegabnds(1),omegabnds(2)); % handle to integrate shape

[shape, shape_int] = generateUncertaintyDist(theta{1}, omegabnds, [1,1]);
ftheta_vec = msubs(shape,theta{1},z3)'; % values of shape over a range of thetas

crispw_func_vec = msubs(w{1},[z{1}],[P*XY])'; % plugs in for z1 and z2 in w. Creates a 2601x1 vector of msspolys that depend on uncertain parameter
crispw_mat = min(1,msubs(crispw_func_vec,theta{1},z3)).^6; % creates a 2601 x 51 matrix that needs to be multiplied by ftheta and summed row-wise to get grid of betas
volume = double(shape_int(msspoly(1)));
beta = (omegabnds(2)-omegabnds(1))*(1/length(z3))*crispw_mat*ftheta_vec/volume; % creates 2601x1 vector of beta values
beta = reshape(beta,size(X)); % creates a grid in XY of betas. scales by max(beta)

%%

% % Shreyas' beta cleanup code
% [S,~] = size(beta) ; % size of beta matrix in rows/cols
% R = S*(v+RT)/(z1bnds(2) - z1bnds(1)); % number of rows/cols for "cleanup
%                      % radius", per unit length
% rrOff = S/(z1bnds(2) - z1bnds(1)) ;
% ccOff = S/(z1bnds(2) - z1bnds(1)) ;
% [rr, cc] = meshgrid(1:S);
% C = sqrt((rr-S/2).^2+(cc-S/2).^2) <= 1.25*R; % TWEAK THE R IN THIS EQUATION
% beta_orig = beta ;
% beta(~C) = 0 ;
%% Then plot it!! Plots beta surface and alpha level-sets
figure(2)
hold on
surf(z1,z2,beta,'EdgeColor','none')

% Shreyas' plot target set on surface
x = linspace(-RT,RT) ;
y = [sqrt(RT.^2-x.^2),-sqrt(RT.^2-x.^2)] ;
plot(x,y(1:end/2),'r',x(end:-1:1),y(end/2+1:end),'r','LineWidth',2)

figure(1)
[C,h] = contour(z1,z2,beta,linspace(0.2,1,5));
tl = clabel(C);

for i = 2:2:length(tl)
   oldLabelText = tl(i).String;
   newLabelText = ['\alpha = ',oldLabelText];
   tl(i).String = newLabelText;
end

%% Plot probability distribution
figure(3)
plot(z3,ftheta_vec/volume,'b')
xlabel('Steering Angle Rate \omega')
ylabel('Probability')
title('\omega Probability Distribution')
axis([omegabnds(1), omegabnds(2), 0, max(ftheta_vec/volume)*1.1])

%% Plot trajectories round 2
TSPAN = [0 T];
Y0 = [-0.5;0]; % starting point

num_traj = 20; % number of trajectories to plot
omegas = zeros(num_traj,1); % initialize omegas
randos = rand([num_traj,1]); % vector of random numbers on (0,1)
probs = cumsum(ftheta_vec*((omegabnds(2) - omegabnds(1))/(length(ftheta_vec)-1)))/volume; % integrates from left omega bound on up

% assign omegas for trajectories through random numbers
for i = 1:length(randos)
    omegas(i) = z3(find(probs>=randos(i),1)) + ((omegabnds(2) - omegabnds(1))/(length(ftheta_vec)-1))/2;
end

figure(1)
hold on

% time to plot dat shit again! 
for i = 1:length(omegas)
    % use ode45 to get trajectories (in X-Y plane)
    [TOUT, YOUT] = ode45(@(t,y) [v*cos(omegas(i)*t); v*sin(omegas(i)*t)], TSPAN, Y0);
    
    plot(YOUT(1,1),YOUT(1,2),'bo') % first point
    plot(YOUT(:,1),YOUT(:,2),'Color',[0.5;0.5;0.5]) % trajectory
    plot(YOUT(end,1),YOUT(end,2),'mx') % last point
end

contour(z1,z2,TargetSet,[0,0],'r'); % plot target set (where hXT is >=0)



