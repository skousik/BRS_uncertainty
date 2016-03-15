% Plotting trajectories for monte carlo

%% Setup
clc
clear
close all

% final time
T = 1 ;
    
% longitudinal velocity
v = 0.5 ;

% starting heading
nominal_theta = 0 ;

% omega bounds (omega lies within +/- omegabd
wBounds = [-pi/4, pi/4] ;
lambda_w = pi/2 ;

% target set radius
RT = sqrt(0.125) ;

% omega distribution
uniform = false ; % whether or not to use uniformly-distributed omegas\
plot_w = true ; % whether or not to plot the omega dist
left_weight = 1 ;
center_weight = 0 ;
right_weight = 1 ;
weights = [left_weight, right_weight] ;

% granularity of X, Y, and omega spaces
Nx = 300 ; % number of x points
Ny = 300 ; % number of y points
Nw = 300 ; % number of omegas

xbds = [-1, 1] ;
ybds = [-1, 1] ;
x = linspace(xbds(1),xbds(2),Nx) ;
y = linspace(ybds(1),ybds(2),Ny) ;

XYbds = [-1, 1 ; -1, 1] ;

%% Determine trajectory endpoints for uniform distribution
w_uniform = linspace(wBounds(1),wBounds(2),Nw) ;
endpoints = findDubinsTrajectoryEndpoints(T, v, w_uniform) ;

%% Find original BRS
% Create grid of points to fill with tallies
mag = zeros(Ny,Nx,Nw) ;

% magical matrix tour
X = repmat(x,Ny,1,Nw) ;
Y = repmat(y',1,Nx,Nw) ;

End1 = repmat(reshape(endpoints(1,:),1,1,Nw),Nx,Ny) ;
End2 = repmat(reshape(endpoints(2,:),1,1,Nw),Nx,Ny) ;

XEnd = X + End1 ;
YEnd = Y + End2 ;

mag = sqrt(XEnd.^2 + YEnd.^2) ;

% get uniform BRS
BRS = sum(mag <= RT, 3) > 0 ;
N_BRS = sum(sum(double(BRS))) ;

%% Create discrete nonuniform distribution of omega
w = msspoly('pa',1) ; % uncertainty
[shape, ~] = generateUncertaintyDist(w, wBounds, weights) ;
scale = lambda_w / (3*Nw) ;
w_distribution = scale*msubs(shape,w, linspace(wBounds(1),wBounds(2),3*Nw)) ;

%% Monte Carlo
% For the w_distribution, find the endpoints of every w
tic
endpoints = findDubinsTrajectoryEndpoints(T, v, w_distribution) ;
toc
%%
endx = endpoints(1,:) ;
endy = endpoints(2,:) ;

%%
% Create an (Nx*Ny) x Nw matrix of indices
idx_mat = zeros(Nx*Ny,Nw) ;

for point = 1:Nx*Ny
    idx = randperm(3*Nw) ;
    idx_mat(point,:) = idx(1:3:end) ;
end

idx_mat = reshape(idx_mat,Ny,Nx,Nw) ;
%%
% create Endx and Endy matrices by using idx_mat to index into endx and
% endy, then reshaping the resulting matrices
Endx = endx(idx_mat) ;
Endy = endy(idx_mat) ;

%%
% Create random-er BRS
X = repmat(x,Ny,1,Nw) ;
Y = repmat(y',1,Nx,Nw) ;

XEnd = X + Endx ;
YEnd = Y + Endy ;

mag = sqrt(XEnd.^2 + YEnd.^2) ;

% get uniform BRS
beta = sum(mag <= RT, 3)./Nw ;

%% Plot beta distribution
figure

subplot(2,2,1)
hold on
surf(x,y,beta,'EdgeColor','none')

% Plot target set
xRT = linspace(-RT,RT) ;
yRT = [sqrt(RT.^2-xRT.^2),-sqrt(RT.^2-xRT.^2)] ;
plot(xRT,yRT(1:end/2),'r',xRT(end:-1:1),yRT(end/2+1:end),'r','LineWidth',2)

title('BRS Colored by \alpha-Level')

axis([-1,1,-1,1])
axis square

%% Plot extra beta distribution to show isometric
subplot(2,2,2)
hold on
surf(x,y,beta,'EdgeColor','none')
view([45 45])

% Plot target set
xRT = linspace(-RT,RT) ;
yRT = [sqrt(RT.^2-xRT.^2),-sqrt(RT.^2-xRT.^2)] ;
plot(xRT,yRT(1:end/2),'r',xRT(end:-1:1),yRT(end/2+1:end),'r','LineWidth',2)

title('BRS Isometric view')

axis([-1,1,-1,1])
axis square

%% Plot omega distribution
if plot_w
    thS = linspace(-omegabd, omegabd, Nw) ;
    scale = 1/double(mu_theta{1}(msspoly(1))) ;
    fplot = ftheta_vec*scale*2*omegabd ;
    
    subplot(2,2,[3,4])
    plot(thS,fplot,'b')
    axis([-omegabd, omegabd, 0, max(fplot)*1.1])
    xlabel('Steering Angle Rate \omega')
    ylabel('Probability')
    title('\omega Probability Distribution')
end
