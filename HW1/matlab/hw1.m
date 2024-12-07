%% Q1
clc; clear;

global_variabels

syms v_an3

r_an_vec = [5000;  8500; 0];
v_an_vec = [2.15; -2.05; v_an3];

a = power(myu * T^2 / (4*pi^2), 1/3);
r_an = norm(r_an_vec);
v_an = sqrt(2*myu*(1/r_an - 1/(2*a)));
v_an_3 = double(sqrt(v_an^2 - v_an_vec(1)^2 - v_an_vec(2)^2));

r_an_vec = [5000;  8500; 0];
v_an_vec = [2.15; -2.05; 5.0124];
h = norm(cross(r_an_vec, v_an_vec));

[a,e,i,small_omega, big_omega,f] = kepler_orbital_elements(r_an_vec,v_an_vec, myu)

%% Q2
clc; clear;

global_variabels

r_an_vec = [5000;  8500; 0];
v_an_vec = [2.15; -2.05; 5.0124];

[a,e,i,small_omega,big_omega,~] = kepler_orbital_elements(r_an_vec,v_an_vec, myu);

rp = a * (1 - e);
vp = sqrt(2 * myu *(1/rp - 1/(2*a)));

rp_vec_peri = [rp; 0; 0];
vp_vec_peri = [0; vp; 0];

rp_vec_ECI = peri2ECI(rp_vec_peri, i, small_omega, big_omega);
vp_vec_ECI = peri2ECI(vp_vec_peri, i, small_omega, big_omega);

state_initial = [rp_vec_ECI; vp_vec_ECI];

time_in_minutes = T;
time_interval = [0 time_in_minutes];

% This is where we integrate the equations of motion.
[t_out, state_out] = ode45(@Satellite, time_interval, state_initial,odeset('RelTol',1e-5));

fig1 = figure ("Name","3D figure of the orbit trajectory along with the Earth drawing",'Position',[100 300 900 500]);
hold all

earth_sphere
plot3(state_out(:,1),state_out(:,2),state_out(:,3), "LineWidth", 2, "Color", "#D95319")

xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
grid on
grid minor
title("3D figure of the orbit trajectory along with the Earth drawing")
subtitle("Almog Dobrescu 214254252")
legend({'The Earth','The orbit'},'FontSize',11 ,'Location','northeast')
%exportgraphics(fig1, 'grap1.png','Resolution',1200);

for i = 1:length(state_out(:,1))
    [a(i),e(i),inclination(i),small_omega(i),big_omega(i),f(i)] = kepler_orbital_elements(state_out(i,1:3),state_out(i,4:6), myu);
end

fig2 = figure ("Name","Orbital Elements",'Position',[300 300 900 500]);
hold all

subplot(5,1,1)

plot(t_out,a, "LineWidth", 1.5)
grid on
ylabel('a [km]')
xlabel('t [sec]')
title("a [km]")
subtitle("Almog Dobrescu 214254252")

subplot(5,1,2)

plot(t_out,e, "LineWidth", 1.5)
grid on
ylabel('e [-]')
xlabel('t [sec]')
title("e [-]")
subtitle("Almog Dobrescu 214254252")

subplot(5,1,3)

plot(t_out,inclination, "LineWidth", 1.5)
grid on
ylabel('i [rad]')
xlabel('t [sec]')
title("i [rad]")
subtitle("Almog Dobrescu 214254252")

subplot(5,1,4)

plot(t_out,small_omega, "LineWidth", 1.5)
grid on
ylabel('ω [rad]')
xlabel('t [sec]')
title("ω [rad]")
subtitle("Almog Dobrescu 214254252")

subplot(5,1,5)

plot(t_out,big_omega, "LineWidth", 1.5)
grid on
ylabel('Ω [rad]')
xlabel('t [sec]')
title("Ω [rad]")
subtitle("Almog Dobrescu 214254252")

%exportgraphics(fig2, 'grap2.png','Resolution',1200);

%% Q3
clc; clear;

global_variabels

r_an_vec = [5000;  8500; 0];
v_an_vec = [2.15; -2.05; 5.0124];

[a,e,i,small_omega,big_omega,~] = kepler_orbital_elements(r_an_vec,v_an_vec, myu);

r = a * (1 - e);
v = sqrt(2 * myu *(1/r - 1/(2*a)));
delta_v = [10e-3;0;0];
f = 0; % perigee
theta = small_omega+f;
h_vector = cross(r_an_vec, v_an_vec);
h = norm(h_vector);
p = h^2/myu;

B2 = [2*a/(r*v)*(2*a-r)                     , 0                           , 0                            ;
      2/v*(e+cos(f))                        , r/(a*v)*sin(f)              , 0                            ;
      0                                     , 0                           , r/h*cos(theta)               ;
      0                                     , 0                           , r/h*sin(theta)/sin(i)        ;
      2/(e*v)*sin(theta)                    , -1/(e*v)*(2*e+r/a*cos(f))   , -r/h*sin(theta)*cos(i)/sin(i);
      -2*sqrt(1-e^2)/(e*v)*(1+e^2*r/p)*sin(f), r*sqrt(1-e^2)/(e*v*a)*cos(f), 0];

B2 * delta_v
 

