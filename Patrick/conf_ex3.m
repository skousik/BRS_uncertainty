% Implementing 1D dynamics for paper
clear all
clc

degree = 14;
T = 1;
freeFinalTime = 0;

thetabnds = 0.5;
xbnds = 1;

t = msspoly('t',1);

x{1} = msspoly('xa',1);
theta{1} = msspoly('pa',1);

xA = x{1};
pA = theta{1};
f{1} = T * [pA;0];

hX{1,1} = [(xbnds-xA(1))*(xA(1)+xbnds)];
hXt{1,1} = [(xbnds-xA(1))*(xA(1)+xbnds);  (theta{1}+thetabnds)*(thetabnds-theta{1})];

hTheta{1,1} = [(theta{1}+thetabnds)*(thetabnds-theta{1})];

hXT{1,1} = [0.25 - (xA-0.5)^2];

dl{1} = boxMoments(x{1},-xbnds,xbnds);
dtheta{1} = boxMoments(theta{1}, -thetabnds, thetabnds);
% % %%
% % mean = 0;
% % sigma = 0.2;
% % offset = 0.025347/2;
% % 
% % A = 1/(sigma*sqrt(2*pi));
% % 
% % syms b
% % shape = A*exp(-((b-mean)^2/(2*sigma^2))) + offset;
% % taylor(shape,'Order',9)
% % %%
% shape = (3509133408902789453125*theta{1}^8)/1729382256910270464 - (140365336356111578125*theta{1}^6)/216172782113783808 + (5614613454244463125*theta{1}^4)/36028797018963968 - (224584538169778525*theta{1}^2)/9007199254740992 + 578589305386791743/288230376151711744

%%%%%% WE THINK THIS SHOULD BE DTHETA IN NEXT LINE
w = solve_BRS_conf(t,x,theta,f,hX,hXT,hTheta,dtheta,dl,degree);

%% Time to plot

% shape = 1;
% shape = -(theta{1}-1)*(theta{1}+1);
% shape = -(theta{1}-.75)*(theta{1}+1.25);
% shape = -(theta{1}-0.5)*(theta{1}+1.5);

shape = -(theta{1}-0.5)*(theta{1}+0.5)
mu_theta{1} = boxMoments_shape(theta{1},shape,-thetabnds,thetabnds);

x1 = linspace(-xbnds,xbnds,81);
x2 = linspace(-thetabnds,thetabnds,51);
[X,Y] = meshgrid(x1,x2);
XY = [X(:) Y(:)]';

subplot(311)
hold on
Wrestr = msubs(w{1},[x{1};theta{1}],XY);
contour(x1,x2,reshape(Wrestr,size(X)),[1,1]);
axis([-xbnds xbnds -0.5 0.5])
% return
subplot(312)
hold on
Wrestr = 1*msubs(mu_theta{1}(w{1}),x{1},x1);
plot(x1,Wrestr/double(mu_theta{1}(msspoly(1))));

% trying to implement beta thing
ftheta_vec = msubs(shape,theta{1},x2);
crispw_vec = min(1,msubs(w{1},[x{1};theta{1}],XY)).^6;
crispw_vec = reshape(crispw_vec,size(X));
beta = (1/length(x2))*ftheta_vec*crispw_vec/double(mu_theta{1}(msspoly(1)));
subplot(313)
hold on
plot(x1,beta);`
