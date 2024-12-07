%% A+B
clc; clear; close all;

global_variabels

T1 = 6e3; % [sec]
T2 = T1;

a1 = (mu*T1^2/(4*pi^2))^(1/3);
a2 = (mu*T2^2/(4*pi^2))^(1/3);

alpha = 0.01*pi/180;
z2_0 = 1;

z2_dot_0_over_n = -sqrt((a1*tan(alpha))^2-z2_0^2);

syms n tau x0 y0

phi = [4-3*cos(n*tau)      , 0, 1/n*sin(n*tau)    , 2/n*(1-cos(n*tau))        ;
       6*(sin(n*tau)-n*tau), 1, 2/n*(cos(n*tau)-1), 1/n*(4*sin(n*tau)-3*n*tau);
       3*n*sin(n*tau)      , 0, cos(n*tau)        , 2*sin(n*tau)              ;
       6*n*(cos(n*tau)-1)  , 0, -2*sin(n*tau)     , 4*cos(n*tau)-3            ];

phi_11 = [phi(1,1), phi(1,2);
          phi(2,1), phi(2,2)];

phi_12 = [phi(1,3), phi(1,4);
          phi(2,3), phi(2,4)];

phi_21 = [phi(3,1), phi(3,2);
          phi(4,1), phi(4,2)];

phi_22 = [phi(3,3), phi(3,4);
          phi(4,3), phi(4,4)];

req = -inv(phi_12)*phi_11*[x0;y0];

final_v = phi_21*[x0;y0]+phi_22*req;

x0 = 0;
y0 = -1;
x0_dot = 0;
y0_dot = 0;
n = sqrt(mu/a1^3);

final_v = simplify(subs(final_v));
eq = final_v(1) == 0;

% t_1 = solve(eq, "ReturnConditions",true);

y_dot_req = -n/6/pi

sqrt(0.74267^2+1^2)

atan2(1,-0.74267)

nt = pi-2.2096

t = nt/n

delta_z = -(-n*sin(0.9320)-0.74267*n*cos(0.9320))

%% C

clc; clear; close all;

global_variabels

T = 6e3; % [sec]
a = (mu*T^2/(4*pi^2))^(1/3);
n = sqrt(mu/a^3);
nt = pi-2.2096;
t_z_eq_0 = nt/n;

x0 = 0;
y0 = -1;
z0 = 1;
x0_dot = 0;
y0_dot = 0;
z0_dot = -0.74267*n;
delta_z = -(-n*sin(0.9320)-0.74267*n*cos(0.9320));
delta_y1 = -5.5556e-5;

ts1 = linspace(0, t_z_eq_0, 1000);
ts2 = linspace(t_z_eq_0, T, 1000);
xs1=[];
ys1=[];
zs1=[];
x_dots1=[];
y_dots1=[];
z_dots1=[];

xs2=[];
ys2=[];
zs2=[];
x_dots2 = [];

for i = 1:length(ts1)
    [xs1(end+1), ys1(end+1), zs1(end+1), x_dots1(end+1), y_dots1(end+1), z_dots1(end+1)] = calculate_CW(x0, y0, z0, x0_dot, y0_dot+delta_y1 , z0_dot, n, ts1(i));
end
for i = 1:length(ts1)
    [xs2(end+1), ys2(end+1), zs2(end+1), x_dots2(end+1), ~, ~] = calculate_CW(xs1(end), ys1(end), zs1(end), x_dots1(end), y_dots1(end), z_dots1(end)+delta_z, n, ts2(i)-ts1(end));
end


fig1 = figure ("Name","3D figure of the orbit trajectory",'Position',[100 300 900 500]);
hold all

plot3(xs1, ys1, zs1, "LineWidth", 2, "Color", "#D95319")
plot3(xs2, ys2, zs2, "LineWidth", 2, "Color", "#D95319", "HandleVisibility","off")
plot3(xs1(1), ys1(1), zs1(1),"hexagram", "LineWidth", 2, "Color", "#0072BD")
plot3(xs1(end), ys1(end), zs1(end),"hexagram", "LineWidth", 2, "Color", "#7E2F8E")
plot3(xs2(end), ys2(end), zs2(end),"hexagram", "LineWidth", 2, "Color", "#FF0000")

xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
grid on
grid minor
title("3D figure of the orbit trajectory along with the Earth drawing")
subtitle("Almog Dobrescu 214254252")
legend({'The orbit','Start', 't = 890[sec]', 'End'},'FontSize',11 ,'Location','northwest')
% exportgraphics(fig1, 'grap1.png','Resolution',600);

%% D
global_variabels

T = 6e3; % [sec]
a = (mu*T^2/(4*pi^2))^(1/3);
n = sqrt(mu/a^3);
delta_x_meas = 1e-3;
delta_y_meas = 1e-3;
delta_z_meas = 1e-3;
delta_x_dot_meas = 1e-5;
delta_y_dot_meas = 1e-5;
delta_z_dot_meas = 1e-5;
ntau = 2*pi;

max_xf_miss = abs(delta_x_meas*(4-3*cos(ntau))) + abs(delta_x_dot_meas/n*sin(ntau)) + abs(2*delta_y_dot_meas/n*(1-cos(ntau)))
max_yf_miss = abs(6*delta_x_meas*(sin(ntau) - ntau)) + abs(delta_y_meas) + abs(2*delta_x_dot_meas/n*(cos(ntau)-1)) + abs(delta_y_dot_meas/n*(4*sin(ntau) -3*ntau))
max_zf_miss = abs(delta_z_meas*cos(ntau)) + abs(delta_z_dot_meas/n*sin(ntau))


