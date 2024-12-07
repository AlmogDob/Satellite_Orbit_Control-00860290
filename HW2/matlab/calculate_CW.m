function [x, y, z, x_dot, y_dot, z_dot] = calculate_CW(x0, y0, z0, x0_dot, y0_dot, z0_dot, n, t)

phi = [4-3*cos(n*t)      , 0, 1/n*sin(n*t)    , 2/n*(1-cos(n*t))        ;
       6*(sin(n*t)-n*t), 1, 2/n*(cos(n*t)-1), 1/n*(4*sin(n*t)-3*n*t);
       3*n*sin(n*t)      , 0, cos(n*t)        , 2*sin(n*t)              ;
       6*n*(cos(n*t)-1)  , 0, -2*sin(n*t)     , 4*cos(n*t)-3            ];

phi_11 = [phi(1,1), phi(1,2);
          phi(2,1), phi(2,2)];
phi_12 = [phi(1,3), phi(1,4);
          phi(2,3), phi(2,4)];
phi_21 = [phi(3,1), phi(3,2);
          phi(4,1), phi(4,2)];
phi_22 = [phi(3,3), phi(3,4);
          phi(4,3), phi(4,4)];
final_pos_xy = phi_11*[x0;y0]+phi_12*[x0_dot;y0_dot];
final_vel_xy = phi_21*[x0;y0]+phi_22*[x0_dot;y0_dot];

x = final_pos_xy(1);
y = final_pos_xy(2);
x_dot = final_vel_xy(1);
y_dot = final_vel_xy(2);

z = z0*cos(n*t) + (z0_dot/n)*sin(n*t);
z_dot = -n*z0*sin(n*t) + z0_dot*cos(n*t);


end

