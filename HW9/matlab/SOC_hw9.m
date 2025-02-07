clc; clear; close all;

%% Part A

mu = 3.986e5; % [km^3/s^2]
R_earth = 6378; % [km] 
T = 23*60*60+56*60+4.09; % [sec]
a = (mu * T^2 / 4 / pi^2)^(1/3); % [km]
n = 2*pi/T; % [1/sec]
f_max = 0.0001e-3; % [km/sec]

lambda_n_deg = 30;
lambda_s_deg = 75;

lambda_n_ddot = -0.00168*sind(2*(lambda_n_deg - lambda_s_deg))

delta_lambda_0 = -sign(lambda_n_ddot)*0.05; % [deg]
delta_lambda = abs(delta_lambda_0);
Delta_L = delta_lambda * pi / 180 * a; % [km]

lambda_0_dot = -2*sign(lambda_n_ddot)*sqrt(-lambda_n_ddot*delta_lambda_0)

tm = 2*abs(lambda_0_dot / lambda_n_ddot)

delta_a = 0.132*sind(2*(lambda_n_deg - lambda_s_deg))*tm

%% Part C
time_interval = 0:1:365*3;

state = [lambda_n_deg; 0; a]; % [lambda, labmbda_dot, a]

 % This is where we integrate the equations of motion.
[t_out, state_out] = ode45(@ODE, time_interval, state, odeset('RelTol',5e-14,'AbsTol',5e-14));

fig1 = figure ("Name","1",'Position',[100 300 900 500]);

plot(t_out/365, state_out(:,1), "LineWidth",1.5)

xlabel('time [year]','FontSize',15,Interpreter='latex')
ylabel('$\lambda$ [deg]','FontSize',15,Interpreter='latex')
grid on
grid minor
title("Drift due to The Tesseral Harmonics")
subtitle("Almog Dobrescu 214254252")
% legend({''},'FontSize',11 ,'Location','northwest')
% exportgraphics(fig1, 'graph1.png','Resolution',300);



%% Functions
function d_state_dt = ODE(t, state)
    mu = 3.986e5; % [km^3/s^2]
    
    lambda = state(1);
    lambda_dot = state(2);
    a = state(3);
    
    n = sqrt(mu/a^3);

    a_dot = 6.06*sind(2*(lambda-27.8))-1.48*sind(4*(lambda-30.28))+0.19*sind(2*(lambda-31.92));

    d_state_dt = [lambda_dot; 3*n/4/a*a_dot; a_dot]; 
end

