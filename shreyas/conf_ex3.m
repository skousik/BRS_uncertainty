% alpha-level confidence
%
% very simple system dynamics (x_dot = uncertain_parameter)

%% Reset
clc
clear
close all
tic
%% User Input
degree = 18 ;
w_crispness = 6 ;
degree_mu_theta_poly = 6 ;
% left_weights = [1] ;
left_weights = [5] ;
w_colors = {[0.9 0.2 0],[0 0.2 0.9]} ;
w_true_colors = {[1 0 0],[0 0 1]} ;
true_theta_granularity = 100 ;

% X bounds
BxaU =  1 ; % lower
BxaL = -1 ; % upper

% mu_theta support bounds
BthU =  0.5 ; % lower
BthL = -0.5 ; % upper
dTh = BthU - BthL ; % Lebesgue size of Theta domain

% Theta overall domain bounds (must be >= than spt(mu_theta))
BthBoxU =  1 ;
BthBoxL = -1 ;
n = 1 ; % degree of theta bounding polynomial

% Target set bounds
BXTL =  0.0 ; % X lower
BXTU =  1.0 ; % X upper
BTTL = -1.0 ; % theta lower
BTTU =  1.0 ; % theta upper

%% Problem Setup
% Create variables:
time  = msspoly('time',1) ;
xa{1} = msspoly('xa',1);
th{1} = msspoly('th',1);

x = xa{1}(1) ;
t = th{1}(1) ;

% Create system dynamics
f{1} = [t; 0] ;

% Define X
hXa{1,1} = [(BxaU-x)*(x-BxaL)] ;

% Support of mu_theta
hTh{1,1} = [(BthU-t)*(t-BthL)];

% Define target set X_T
hXT{1,1} = [-1*(x-BXTU)*(x-BXTL);
            -1*(t-BTTU)*(t-BTTL)];
        
% Integrator for w
dl{1} = boxMoments(xa{1}, BxaL, BxaU) ;

% Theta set
u = BthBoxU ;
l = BthBoxL ;
dth{1} = boxMoments(th{1}, l, u) ;

%% Optimize!
tic
[w,sol] = solve_BRS_conf(time, xa, th, f, hXa, hXT, hTh, dth, dl, degree) ;
toc
%% Plottin'
tic
close all
xaS = linspace(BxaL, BxaU, 100) ; % x space
thS = linspace(BthBoxL, BthBoxU, 100) ; % theta space
[Xa, Th] = meshgrid(xaS, thS) ;
XaTh = [Xa(:) Th(:)]';
orange = [0.8 0.5 0] ;
outer_green = [0.1 0.5 0] ;
green = [0.1 0.7 0] ;

% subplot(311)
figure(1)
grid on
hold on

% Plot X_T
plot([0,1,1,0,0],2.*[-0.5,-0.5,0.5,0.5,-0.5],'Color',orange,'LineWidth',1.5)

% Plot the true BRS
xa_true = [-0.5, 0.5, 1.0,  1.0,  0.5, -0.5] ;
th_true = [ 0.5, 0.5, 0.0, -0.5, -0.5,  0.5] ;
plot(xa_true, th_true, 'Color', outer_green, 'LineWidth', 2) ;

% Plot the estimate
W1 = msubs(w{1},[xa{1};th{1}],XaTh);
contour(xaS,thS,reshape(W1,size(Xa)),[1,0], ...
        'Color', green, 'LineStyle','--', 'LineWidth', 2);

xlabel('x')
ylabel('\theta')

set(gca,'FontSize',20)

%% Plotting mu_theta and beta
% Here, we approximate mu_theta with an even-ordered polynomial on its
% support, and vary the distribution from left- to right-tailed by changing
% the weights of the repeated zeros at the left and right bounds of
% spt(mu_theta)

% Set the variables that don't change through loops
gran = true_theta_granularity ;
xaS = linspace(BxaL, BxaU, gran) ; % x space
thS = linspace(BthL, BthU, gran) ; % theta space
[Xa, Th] = meshgrid(xaS, thS) ;
XaTh = [Xa(:) Th(:)]';

% Create matrix mesh representing w
O = ones(gran/2) ;
Z = zeros(gran/2) ;
W = [triu(O),  O, triu(O)' ; Z, triu(O), O] ;
removeCols = 3:3:1.5*gran ;
W(:,removeCols) = [] ;

% colors for plots
C = w_colors ;
C_true = w_true_colors ;

for i = 1:length(left_weights)
  % Shape of mu_theta
    % weight of left root
    l_wt  = left_weights(i) ;
    % weight of right root
    r_wt  = degree_mu_theta_poly - l_wt ;
    % RN-derivative shape function as a polynomial
    [shape, ~] = generateUncertaintyDist(th, ...
                            [BthL, BthU], [l_wt, r_wt]) ;

    mu{1} = boxMoments_shape(th{1}, shape, BthL, BthU) ;

  % Create mu_theta distribution and integral of w across Theta
    k = w_crispness ; % determines crispness of w plot
    
    % f_muth is the RN-derivative
    f_muth_1 = msubs(shape,th{1},thS) ;
    
    % w_xath is w as a discretized function of x, integrated over Theta
    w_xath = reshape(min(1,(msubs(w{1},[xa{1}; th{1}], XaTh))).^(2^k), size(Xa)) ;
    
    % beta is the discretized w evaluated at each x
    beta = (dTh/length(thS))*f_muth_1*w_xath ;
    
    
    [shape, ~] = generateUncertaintyDist(th, ...
                        [BthL, BthU], [r_wt, l_wt]) ;
                        % f_muth is the RN-derivative
    f_muth_2 = msubs(shape,th{1},thS) ;
  
  % Generate true beta
    F = repmat(f_muth_2',1,gran)*(dTh/length(thS)) ;
    intW = sum(F.*W) ;
    
  % Plot w
%     subplot(312)
    figure(2)
    hold on
    grid on
    plot(xaS,beta,'--','Color',C{i},'LineWidth',2);
    hold on
    xlabel('x')
    ylabel('Probability of success')
    set(gca,'FontSize',20)
    
  % Plot true w
    plot(linspace(-0.5,1,gran),intW,'Color',C_true{i},'LineWidth',2);
    
  % Plot mu_theta distribution
%     subplot(313)
    figure(3)
    hold on
    plot(thS,f_muth_1,'Color',C{i},'LineWidth',2);
    grid on
    xlabel('\theta')
    ylabel('\mu_\theta')
    set(gca,'FontSize',24)
end
toc
% %% Single distribution, calculate true vs. estimate
% l_wt = 1 ;
% r_wt = degree_mu_theta_poly - l_wt ;
% shape = -1*(t-BthU)^(r_wt)*(t-BthL)^(l_wt) ; % f as a function of theta, with weights flipped
% gran = 600 ; % fineness of true w plot
% thetas = linspace(BthL,BthU,gran) ; % range of thetas
% f = msubs(shape,t, thetas) ; % shape discretized across Theta
% f = f'/sum(f) ; % normalized, discretized shape over Theta
% F = repmat(f,1,gran) ; % matrix of repeated f
% 
% intW = sum(F.*W) ;
% x_vec = linspace(-0.5,1,gran) ;
% 
% % Plot true curve
% subplot(224)
% hold on
% grid on
% plot(x_vec, intW, 'Color', orange, 'LineWidth', 1.5) ;
% 
% % Create mu_theta distribution and integral of w across Theta
% mu{1} = boxMoments_shape(th{1}, shape, BthL, BthU) ;
% scale = 1/double(mu{1}(msspoly(1))) ;
% k = w_crispness ; % determines crispness of w plot
% 
% % Reset shape in order to create the right version
% shape = -1*(t-BthU)^(l_wt)*(t-BthL)^(r_wt) ; % f as a function of theta
% 
% % f_muth is the RN-derivative
% f_muth = msubs(shape,th{1},thS) ;
% 
% % w_xath is w as a discretized function of x, integrated over Theta
% w_xath = reshape(min(1,(msubs(w{1},[xa{1}; th{1}], XaTh))).^(2^k), size(Xa)) ;
% max(w_xath)
% 
% % beta is the discretized w evaluated at each x
% beta = (1/length(thS))*f_muth*w_xath*scale ;
% beta = beta./max(beta) ;
% 
% % Create gamma(x)
% x_vec = [ linspace(-0.5,0.5,66), linspace(0.5,1,34)] ;
% t_vec = [-linspace(-0.5,0.5,66), 1-linspace(0.5,1,34)] ;
% 
% % Plot w
% subplot(224)
% plot(xaS,beta,'--','Color',[0,0,0.5],'LineWidth',1.5);
% hold on
% xlabel('X')
% ylabel('\intw(x,\theta)d\mu_{\theta}')
% title('\alpha-level as a function of x \in X')
% legend('True','Estimate','Location','northwest')
% 
% % Plot mu_theta distribution
% h = subplot(221) ;
% hold on
% plot(f_muth*scale,thS,'b')
% grid on
% p = get(h,'position') ;
% p(3) = p(3)/2 ;
% p(1) = p(1)*2.2 ;
% set(h, 'position', p) ;
% ylabel('\Theta')
% xlabel('\mu_\theta')
% title('Distribution of \mu_\theta')
% axis([0,3,-1,1])



% %% Original w plot
% thS = linspace(BthL, BthU, 100) ; % theta space
% subplot(312)
% hold on
% scale = 1/double(mu{1}(msspoly(1)))
% Wrestr = msubs(mu{1}(w{1}),xa{1},xaS)*scale;
% plot(xaS,Wrestr);
% xlabel('X')
% ylabel('\intw(x,\theta)d\mu_{\theta}')


