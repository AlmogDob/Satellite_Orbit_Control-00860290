mu = 3.986e5; % [km^3/s^2]
R_earth = 6378; % [km] 
T = 6e3; %[sec]
e = 0;
a = (mu*T^2/4/pi^2)^(1/3); %[km]
n = sqrt(mu/a^3);
t_f = 2000; %[sec]
accel_limit = 8e-6; %[km/sec^2]

pos0_vec = [0;-1;1]; %[km]
v0_vec = [0;0;-0.74267*n]; %[km/sec]
posf_vec = [0;0;0];
vf_vec = [0;0;0];

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

syms tau
phi = [4-3*cos(n*tau)      , 0, 0            , 1/n*sin(n*tau)    , 2/n*(1-cos(n*tau))        , 0             ;
       6*(sin(n*tau)-n*tau), 1, 0            , 2/n*(cos(n*tau)-1), 1/n*(4*sin(n*tau)-3*n*tau), 0             ;
       0                   , 0, cos(n*tau)   , 0                 , 0                         , 1/n*sin(n*tau);
       3*n*sin(n*tau)      , 0, 0            , cos(n*tau)        , 2*sin(n*tau)              , 0             ;
       6*n*(cos(n*tau)-1)  , 0, 0            , -2*sin(n*tau)     , 4*cos(n*tau)-3            , 0             ;
       0                   , 0, -n*sin(n*tau), 0                 , 0                         , cos(n*tau)    ];
phi_11 = [phi(1,1), phi(1,2), phi(1,3);
          phi(2,1), phi(2,2), phi(2,3);
          phi(3,1), phi(3,2), phi(3,3)];
phi_12 = [phi(1,4), phi(1,5), phi(1,6);
          phi(2,4), phi(2,5), phi(2,6);
          phi(3,4), phi(3,5), phi(3,6)];
phi_21 = [phi(4,1), phi(4,2), phi(4,3);
          phi(5,1), phi(5,2), phi(5,3);
          phi(6,1), phi(6,2), phi(6,3)];
phi_22 = [phi(4,4), phi(4,5), phi(4,6);
          phi(5,4), phi(5,5), phi(5,6);
          phi(6,4), phi(6,5), phi(6,6)];

C_star = -inv(phi_12)*phi_11;