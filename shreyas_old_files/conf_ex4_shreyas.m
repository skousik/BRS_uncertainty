% Implementing Dubin's car model for paper
clear all
close all
clc
degree = 9;
T = 1;
freeFinalTime = 0;

nominal_omega = 0;
omegabnds = [nominal_omega - 0.1;nominal_omega + 0.1];
initial_heading = 0 ;
v = 0.5; % speed in m/s
RT = sqrt(0.125) ; % radius of target set
z1Off = 0.0 ;
z2Off = 0.0 ;


z1bnds = [-1;1] ;
z2bnds = [-1;1] ;

% define optimization variables
t = msspoly('t',1);

z{1} = msspoly('za',2);
theta{1} = msspoly('pa',1);

zA = z{1};
pA = theta{1};
f{1} = T*[v - zA(2)*pA; zA(1)*pA; 0];

hX{1,1} = [-(zA(1) - z1bnds(1))*(zA(1) - z1bnds(2)); -(zA(2) - z2bnds(1))*(zA(2)-z2bnds(2))];

hTheta{1,1} = [-(theta{1} - omegabnds(1))*(theta{1} - omegabnds(2))];

hXT{1,1} = [RT^2 - (zA(1)-z1Off)^2 - (zA(2)-z2Off)^2];

dl{1} = boxMoments(z{1}, [z1bnds(1);z2bnds(1)], [z1bnds(2);z2bnds(2)]);
dtheta{1} = boxMoments(pA, omegabnds(1), omegabnds(2));

% optimize
[w,sol] = solve_BRS_conf(t,z,theta,f,hX,hXT,hTheta,dtheta,dl,degree);

%% Discretizing space

% z1 = linspace(z1bnds(1),z1bnds(2),21);
z1 = linspace(z1bnds(1),z1bnds(2),51);
z2 = linspace(z2bnds(1),z2bnds(2),51);
z3 = linspace(omegabnds(1),omegabnds(2),51);

% z4 = linspace(-thetabnds,thetabnds,11);

[X,Y] = meshgrid(z1,z2);
XY = [X(:) Y(:)]';

P = [cosd(initial_heading)  sind(initial_heading);
     sind(initial_heading) -cosd(initial_heading)];

P2 = @(initial_heading) [cosd(initial_heading)  sind(initial_heading);
                        sind(initial_heading) -cosd(initial_heading)];

 
% figure(1)
% hold on
Wrestr = msubs(w{1},[z{1};theta{1}],[XY; nominal_omega*ones(size(X(:)))']);
Wrestr_CARTESIAN = msubs(w{1},[z{1};theta{1}],[P*XY; nominal_omega*ones(size(X(:)))']);
Wrestr = reshape(Wrestr,size(X));
Wrestr_CARTESIAN = reshape(Wrestr_CARTESIAN,size(X));

% contour(z1,z2,Wrestr,[1,1],'b')
% contour(z1,z2,Wrestr_CARTESIAN,[1,1],'b');
TargetSet = msubs(hXT{1,1},z{1},P*XY);
TargetSet = reshape(TargetSet,size(X));
% contour(z1,z2,TargetSet,[0,0],'r');
% axis square

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
l_wt = 1 ;
r_wt = 1 ;

shape = -(theta{1}-omegabnds(1)).^(r_wt)*(theta{1}-omegabnds(2)).^(l_wt);
mu_theta{1} = boxMoments_shape(theta{1},shape,omegabnds(1),omegabnds(2));
% Z = repmat(z3,length(X),1);
% Zstar = [Z(:)]';


ftheta_vec = msubs(shape,theta{1},z3)';
% Zstar = repmat(ftheta_vec,1,length(X));
% Zstar = Zstar(:)'
% ftheta_vec = ones(size(X));
crispw_func_vec = msubs(w{1},[z{1}],[P*XY])';
crispw_mat = min(1,msubs(crispw_func_vec,theta{1},z3).^6);
beta = (1/length(z3))*crispw_mat*ftheta_vec/double(mu_theta{1}(msspoly(1)));
beta = reshape(beta,size(X))    ;

%% beta cleanup
[S,~] = size(beta) ; % size of beta matrix in rows/cols
R = S*(v+RT)/(z1bnds(2) - z1bnds(1)); % number of rows/cols for "cleanup
                                         % radius", per unit length
rrOff = S*z1Off/(z1bnds(2) - z1bnds(1)) ;
ccOff = S*z2Off/(z1bnds(2) - z1bnds(1)) ;

[rr, cc] = meshgrid(1:S);
C = sqrt((rr-S/2).^2+(cc-S/2).^2) <= R; % TWEAK THE R IN THIS EQUATION
beta_orig = beta ;
beta(~C) = 0 ;
% surf(z1,z2,double(C))
subplot(1,3,1)
surf(z1,z2,beta)
subplot(1,3,2)
surf(z1,z2,double(~C))
subplot(1,3,3)
surf(z1,z2,beta_orig)
%%

% % beta = (1/length(ftheta_vec))*ftheta_vec*crispw_vec/double(mu_theta{1}(msspoly(1)));
% beta = (1/length(ftheta_vec))*ftheta_vec*crispw_vec;
% 
% figure(2)
% surf(z1,z2,beta)

%%
% figure(3)
% contour(z1,z2,TargetSet,[0,0],'r');
% hold on
% [C, h] = contour(z1,z2,beta,linspace(0,1.25,6),'b') ;
% clabel(C,h)
% axis square
