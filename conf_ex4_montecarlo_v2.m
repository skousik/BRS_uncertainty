% final time
T = 1 ;

% longitudinal velocity
v = 0.5 ;

% starting heading
nominal_theta = 0 ;

% omega bounds (omega lies within +/- omegabd)
omegabd = 3*pi/4 ;

% target set radius
RT = sqrt(0.125) ;

% omega distribution
uniform = false ; % whether or not to use uniformly-distributed omegas\
plot_omega = true ; % whether or not to plot the omega dist
left_weight = 1 ;
center_weight = 0 ;
right_weight = 1 ;

% granularity of X, Y, and omega spaces
Nx = 300 ; % number of x points
Ny = 300 ; % number of y points
Nw = 100 ; % number of omegas

xbds = [-1, 1] ;
ybds = [-1, 1] ;
x = linspace(xbds(1),xbds(2),Nx) ;
y = linspace(ybds(1),ybds(2),Ny) ;

XYbds = [-1, 1 ; -1, 1] ;


%% BRS Finding
% First, find the binary BRS
omega_distribution = linspace(-omegabd, omegabd, Nw) ;

%% Create distribution of omega
% start with uniform distribution
omega_distribution = linspace(-omegabd, omegabd, Nw) ;

if ~uniform
    % create mss poly to evaluate this thing
    theta{1} = msspoly('pa',1) ; % uncertainty

    l_wt = left_weight ;
    c_wt = center_weight ;
    r_wt = right_weight ;
     % shape determines distribution of uncertainty
    shape = -(theta{1})^(c_wt)*(theta{1}-omegabd)^(l_wt)*(theta{1}+omegabd)^(r_wt);

     % handle to integrate shape
    mu_theta{1} = boxMoments_shape(theta{1},shape,-omegabd,omegabd);

     % values of shape over a range of thetas
    ftheta_vec = msubs(shape,theta{1},omega_distribution)';

    volume = double(mu_theta{1}(msspoly(1)));

    num_traj = Nw;
    omegas = zeros(num_traj,1);
    % 
    randos = rand([num_traj,1]);
    probs = cumsum(ftheta_vec*((2*omegabd)/(length(ftheta_vec)-1)))/volume;

    for i = 1:length(randos)
        omegas(i) = omega_distribution(find(probs>=randos(i),1)) + ...
                                       ((2*omegabd)/(length(ftheta_vec)-1))/2;
    end

    omega_distribution = omegas ;
end

%% Determine trajectory endpoints


%% Monte Carlo sim
% Create grid of points to fill with tallies
mag = zeros(Ny,Nx,Nw) ;

%% magical matrix tour
Z1 = repmat(x,Ny,1,Nw) ;
Z2 = repmat(y',1,Nx,Nw) ;

End1 = repmat(reshape(endpoints(1,:),1,1,Nw),Nx,Ny) ;
End2 = repmat(reshape(endpoints(2,:),1,1,Nw),Nx,Ny) ;

Z1End = Z1 + End1 ;
Z2End = Z2 + End2 ;

mag = sqrt(Z1End.^2 + Z2End.^2) ;

beta = sum(mag <= RT, 3) ;
beta = beta/Nw ;
toc
%% Plot beta distribution
tic % plotting timer
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
if plot_omega
    thS = linspace(-omegabd, omegabd, Nw) ;
    scale = 1/double(mu_theta{1}(msspoly(1))) ;
    fplot = ftheta_vec/sum(ftheta_vec) ;
    
    subplot(2,2,[3,4])
    plot(thS,ftheta_vec/sum(ftheta_vec),'b')
    axis([-omegabd, omegabd, 0, max(fplot)*1.1])
    xlabel('Steering Angle Rate \omega')
    ylabel('Probability')
    title('\omega Probability Distribution')
end
toc