%% Implementing Dubin's car model for paper
clear all
clc
degree = 12;
T = 1;
freeFinalTime = 0;

nominal_omega = 0; % where the distribution is centered around
omegabnds = [nominal_omega - 3*pi/4;nominal_omega + 3*pi/4]; % bounds of distribution
v = 0.5; % speed in m/s

z1bnds = [-1;1]; % bounds on domain of z1
z2bnds = [-1;1]; % bounds on domain of z2

%% define optimization variables
t = msspoly('t',1); % time

z{1} = msspoly('za',2); % states z1 and z2
theta{1} = msspoly('pa',1); % uncertainty

zA = z{1};
pA = theta{1};
f{1} = T*[v - zA(2)*pA;
              zA(1)*pA;
                     0]; % dynamics

hX{1,1} = [-(zA(1) - z1bnds(1))*(zA(1) - z1bnds(2)); -(zA(2) - z2bnds(1))*(zA(2)-z2bnds(2))]; % quadratics define domain of z

hTheta{1,1} = [-(theta{1} - omegabnds(1))*(theta{1} - omegabnds(2))]; % domain of uncertainty

hXT{1,1} = [0.125 - zA'*zA]; % target set is a ball around the origin

dl{1} = boxMoments(z{1}, [z1bnds(1);z2bnds(1)], [z1bnds(2);z2bnds(2)]);
dtheta{1} = boxMoments(pA, omegabnds(1), omegabnds(2));

% optimize
[w,sol] = solve_BRS_conf(t,z,theta,f,hX,hXT,hTheta,dtheta,dl,degree);

%% Discretizing space

z1 = linspace(z1bnds(1),z1bnds(2),101); % discretize z1
z2 = linspace(z2bnds(1),z2bnds(2),101); % discretize z2
z3 = linspace(omegabnds(1),omegabnds(2),101); % discretize z3

[X,Y] = meshgrid(z1,z2); % create a grid of z1 and z2
XY = [X(:) Y(:)]'; % this turns the grid into a long vector
P = [1 0; 0 -1]; % transformation matrix from Cartesian to z-coordinates at steering angle = 0 (corresponds to time = 0)
figure(1)
hold on
Wrestr = msubs(w{1},[z{1};theta{1}],[XY; nominal_omega*ones(size(X(:)))']); % assume uncertainty is equal to nominal uncertainty.
Wrestr_CARTESIAN = msubs(w{1},[z{1};theta{1}],[P*XY; nominal_omega*ones(size(X(:)))']); % this gives you values of w for a grid in X-Y plane
Wrestr = reshape(Wrestr,size(X)); % reshape vector into grid in Z
Wrestr_CARTESIAN = reshape(Wrestr_CARTESIAN,size(X)); % reshape vector into grid in X-Y plane

% contour(z1,z2,Wrestr,[1,1],'b')
% contour(z1,z2,Wrestr_CARTESIAN,[1,1],'b'); % plot Cartesian BRS
TargetSet = msubs(hXT{1,1},z{1},XY); % define target set over grid
TargetSet = reshape(TargetSet,size(X));
contour(z1,z2,TargetSet,[0,0],'r'); % plot target set (where hXT is >=0)

% %% Plot trajectories
% TSPAN = [0 T];
% Y0 = [-1; 0];
% omega = nominal_theta;
% 
% N = 6; %generate <= N^2 points within ROA
% 
% % Get point within ROA
% 
% [z2_ROA, z1_ROA] = find(Wrestr>=1);
% min_z1 = z1(min(z1_ROA(3:end-11)));
% max_z1 = z1(max(z1_ROA(11:end-2)));
% min_z2 = z2(min(z2_ROA(11:end-11)));
% max_z2 = z2(max(z2_ROA(11:end-11)));
% 
% z1plot = linspace(min_z1,max_z1,N);
% z2plot = linspace(min_z2,max_z2,N);
% hold on
% 
% % Z1 = repmat(z1plot,N,1);
% % Z2 = repmat(z2plot,N,1);
% % new_Wrestr = msubs(w{1},[z{1};theta{1}], [Z1(:) Z2(:) nominal_theta*ones(size(Z1(:)))]');
% % z1vals = Z1(new_Wrestr >= 1);
% % z2vals = Z2(new_Wrestr >= 1);
% 
% for i = 1:length(z1plot)
%     for j = 1:length(z2plot)
%     Y0 = [z1plot(i); z2plot(j)];
%     [TOUT, YOUT] = ode45(@(t,y) [v - y(2)*omega; y(1)*omega], TSPAN, Y0);
%     
%     heading = TOUT*nominal_theta;
%     for k = 1:length(TOUT)
%         Q = [cos(heading(k)) sin(heading(k)); sin(heading(k)) -cos(heading(k))];
%         YOUT(k,:) = (Q\YOUT(k,:)')';
%     end
%     
%     plot(YOUT(1,1),YOUT(1,2),'bo')
%     plot(YOUT(:,1),YOUT(:,2),'Color',[0.5;0.5;0.5])
%     plot(YOUT(end,1),YOUT(end,2),'mx')
%     end
% end

%% Time to integrate dat shit!!
l_wt = 1;
r_wt = 1;
shape = -(theta{1}-omegabnds(1))^r_wt*(theta{1}-omegabnds(2))^l_wt; % shape determines distribution of uncertainty
mu_theta{1} = boxMoments_shape(theta{1},shape,omegabnds(1),omegabnds(2)); % handle to integrate shape

ftheta_vec = msubs(shape,theta{1},z3)'; % values of shape over a range of thetas

crispw_func_vec = msubs(w{1},[z{1}],[P*XY])'; % plugs in for z1 and z2 in w. Creates a 2601x1 vector of msspolys that depend on uncertain parameter
crispw_mat = min(1,msubs(crispw_func_vec,theta{1},z3).^6); % creates a 2601 x 51 matrix that needs to be multiplied by ftheta and summed row-wise to get grid of betas
volume = double(mu_theta{1}(msspoly(1)));
beta = (1/length(z3))*crispw_mat*ftheta_vec/volume; % creates 2601x1 vector of beta values
beta = reshape(beta,size(X))./max(beta); % creates a grid in XY of betas. scales by max(beta)

figure(2)
surf(z1,z2,beta,'EdgeColor','none')
figure(1)
[C,h] = contour(z1,z2,beta,linspace(0.2,1,5));
tl = clabel(C);

for i = 2:2:length(tl)
   oldLabelText = tl(i).String;
   newLabelText = ['\alpha = ',oldLabelText];
   tl(i).String = newLabelText
end
%% Plot trajectories round 2
TSPAN = [0 T];
Y0 = [-0.75;-0.25];
% omegas = linspace(omegabnds(1),omegabnds(2),5);
num_traj = 20;
omegas = zeros(num_traj,1);
% 
randos = rand([num_traj,1]);
probs = cumsum(ftheta_vec*((omegabnds(2) - omegabnds(1))/(length(ftheta_vec)-1)))/volume;

for i = 1:length(randos)
    omegas(i) = z3(find(probs>=randos(i),1)) + ((omegabnds(2) - omegabnds(1))/(length(ftheta_vec)-1))/2;
end
% omegas = 
% 
% N = 6; %generate <= N^2 points within ROA
% 
% % Get point within ROA
% 
% [z2_ROA, z1_ROA] = find(Wrestr>=1);
% min_z1 = z1(min(z1_ROA(3:end-11)));
% max_z1 = z1(max(z1_ROA(11:end-2)));
% min_z2 = z2(min(z2_ROA(11:end-11)));
% max_z2 = z2(max(z2_ROA(11:end-11)));
% 
% z1plot = linspace(min_z1,max_z1,N);
% z2plot = linspace(min_z2,max_z2,N);
figure(1)
hold on
% 
% % Z1 = repmat(z1plot,N,1);
% % Z2 = repmat(z2plot,N,1);
% % new_Wrestr = msubs(w{1},[z{1};theta{1}], [Z1(:) Z2(:) nominal_theta*ones(size(Z1(:)))]');
% % z1vals = Z1(new_Wrestr >= 1);
% % z2vals = Z2(new_Wrestr >= 1);
% 
for i = 1:length(omegas)
    [TOUT, YOUT] = ode45(@(t,y) [v*cos(omegas(i)*t); v*sin(omegas(i)*t)], TSPAN, Y0);
    
    plot(YOUT(1,1),YOUT(1,2),'bo')
    plot(YOUT(:,1),YOUT(:,2),'Color',[0.5;0.5;0.5])
    plot(YOUT(end,1),YOUT(end,2),'mx')
end




