mu = 3.986e5; % [km^3/s^2]
R_earth = 6378; % [km] 
T = 6e3; %[sec]
e = 0;
a = (mu*T^2/4/pi^2)^(1/3); %[km]
n = 2*pi/T;
f_max = 4e-5;
pos0_vec = [0;-1;1]; %[km]
v0_vec = [0;0;-7.7765e-4]; %[km/sec]
posf_vec = [0;0;0];
vf_vec = [0;0;0];
x0_vec = [pos0_vec(1);v0_vec(1);pos0_vec(2);v0_vec(2);pos0_vec(3);v0_vec(3)];
xf_vec = [posf_vec(1);vf_vec(1);posf_vec(2);vf_vec(2);posf_vec(3);vf_vec(3)];
t_f = 2000; % [sec]

F = [0    , 1   , 0, 0  , 0   , 0;
     3*n^2, 0   , 0, 2*n, 0   , 0;
     0    , 0   , 0, 1  , 0   , 0;
     0    , -2*n, 0, 0  , 0   , 0;
     0    , 0   , 0, 0  , 0   , 1;
     0    , 0   , 0, 0  , -n^2, 0];
G = [0, 0, 0;
     1, 0, 0;
     0, 0, 0;
     0, 1, 0;
     0, 0, 0;
     0, 0, 1];
Q = [1/(1e-3)^2, 0         , 0         , 0         , 0         , 0         ;
     0         , 1/(1e-5)^2, 0         , 0         , 0         , 0         ;
     0         , 0         , 1/(1e-3)^2, 0         , 0         , 0         ;
     0         , 0         , 0         , 1/(1e-5)^2, 0         , 0         ;
     0         , 0         , 0         , 0         , 1/(1e-3)^2, 0         ; 
     0         , 0         , 0         , 0         , 0         , 1/(1e-5)^2];
R = 1/f_max^2;
Psi = eye(6);


factor = 600;
p = 10*[-n+1i*n/factor; -n-n*1i/factor; -4*n+3*n*1i/factor; -4*n-3*n*1i/factor; -3*n+n*1i/factor; -3*n-n*1i/factor];
% p1=[-1/200+1i*n,-1/200-1i*n];
% p=[p1,4*p1,p1];
K4 = place(F, G, p);

K6 = lqr(F, G, Q, R);

