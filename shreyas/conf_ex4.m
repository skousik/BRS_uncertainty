%% dubins car implementation
close all
clear all
clc

%% User Input
degree = 10 ;
v = 1 ; % speed in m/s
th_nom = 45 ; % heading angle in deg
R = 0.5 ; % radius of target set

% Theta bounds
BthL = -0.01 ; 
BthU =  0.01 ;

% Z Bounds
Bz1L = -1 ;
Bz1U =  1 ;
Bz2L = -1 ;
Bz2U =  1 ;

%% Problem setup
% Variables
time  = msspoly('time',1) ;
za{1} = msspoly('za',2);
th{1} = msspoly('th',1);

z = za{1} ;
t = th{1} ;
th_nom = degtorad(th_nom) ;

% System dynamics
f{1} = [v - z(2)*t ;
        z(1)*t ;
        0] ;
    
% Define Z
hZ{1,1} = [-(z(1)-Bz1L)*(z(1)-Bz1U);
           -(z(2)-Bz2L)*(z(2)-Bz2U)] ;

% Define Theta
hTh{1,1} = [-(t-BthL)*(t-BthU)];

% Define target set
hZT{1,1} = [R - z'*z] ;

% l for cost function
dl{1} = boxMoments(z, [Bz1L;Bz2L], [Bz2L;Bz2U]) ;
dtheta{1} = boxMoments(t, BthL, BthU) ;

%% Optimize!!
w = solve_BRS_conf(time,za,th,f,hZ,hZT,hTh,dtheta,dl,degree) ;

%% Plot BRS
z1 = linspace(Bz1L,Bz1U,51);
z2 = linspace(Bz2L,Bz2U,51);

[Z1, Z2] = meshgrid(z1,z2) ;
Z1Z2 = [Z1(:) Z2(:)]' ;

% Transformation matrix from Z to Cartesian
P = [cos(th_nom)  sin(th_nom);
     sin(th_nom) -sin(th_nom)];
 
% w in z domain
W_Z = msubs(w{1}, [za{1}; th{1}], [Z1Z2; th_nom*ones(size(Z1(:)))']) ;
W_Z = reshape(W_Z, size(Z1));

% w in Cartesian coords
W_C = msubs(w{1}, [za{1}; th{1}], [P*Z1Z2; th_nom*ones(size(Z1(:)))']) ;
W_C = reshape(W_C, size(Z1)) ;

% Plot BRS
contour(z1,z2,W_C,[1,1],'b') ;

%% Plot target set
Tset = msubs(hZT{1,1}, z, Z1Z2) ;
Tset = reshape(Tset, size(Z1)) ;
contour(z1,z2,Tset, [0,0], 'r') ;

axis square