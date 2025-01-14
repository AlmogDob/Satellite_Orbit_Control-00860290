%% A
clc; clear; close all;

global_variabels

pos0_vec = [0;-1;1]; %[km]
v0_vec = [0;0;-0.74267*n]; %[km/sec]
posf_vec = [0;0;0];
vf_vec = [0;0;0];

x0_vec = [pos0_vec(1);v0_vec(1);pos0_vec(2);v0_vec(2);pos0_vec(3);v0_vec(3)];
xf_vec = [posf_vec(1);vf_vec(1);posf_vec(2);vf_vec(2);posf_vec(3);vf_vec(3)];
t_f = 2000;

v_ref = 3e-5; % [m/sec]
rpos0_vec = [0;0;v_ref*t_f];
rv0_vec = -rpos0_vec/norm(rpos0_vec)*v_ref;

xr0_vec = [rpos0_vec(1);rv0_vec(1);rpos0_vec(2);rv0_vec(2);rpos0_vec(3);rv0_vec(3)];

state = [x0_vec-xr0_vec;xr0_vec];
time_interval = [0:0.1:t_f];

 % This is where we integrate the equations of motion.
[t_out, state_out] = ode45(@Satellite, time_interval, state, odeset('RelTol',5e-14,'AbsTol',5e-14));

miss_distance = norm([state_out(end,1),state_out(end,2),state_out(end,3)]);
miss_velocity = norm([state_out(end,4),state_out(end,5),state_out(end,6)]);
x_vec_out = [state_out(:,1) + state_out(:,7),state_out(:,3) + state_out(:,9),state_out(:,5) + state_out(:,11)];
vec_delx_out = [state_out(:,1),state_out(:,3),state_out(:,5)];
x_ref_vec_out = [state_out(:,7),state_out(:,9),state_out(:,11)];

fig1 = figure ("Name","3D Figure of The Orbit Trajectory",'Position',[100 300 900 500]);

plot3(x_vec_out(:,1),x_vec_out(:,2),x_vec_out(:,3), "LineWidth", 2, "Color", cool(1)/1.5)
hold all
plot3(state_out(:,7),state_out(:,9),state_out(:,11), "LineWidth", 2, "Color", "r")
plot3(x_vec_out(1,1),x_vec_out(1,2),x_vec_out(1,3), "o", "LineWidth", 4, "Color", "k")
plot3(x_vec_out(end,1),x_vec_out(end,2),x_vec_out(end,3),'o', "LineWidth", 1, "Color", "m")

xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
grid on
grid minor
title(sprintf("3D Figure of The Orbit Trajectory | miss distance = %g[m] | miss velocity = %g[cm/sec]", miss_distance*10^3, miss_velocity*10^5))
subtitle("Almog Dobrescu 214254252")
legend({'the orbit', 'target trajectory', 'init position', 'end position'},'FontSize',11 ,'Location','northwest')
% exportgraphics(fig1, 'graph1.png','Resolution',600);
%%
colors = cool(9);
%finding the time when norm(vec_delx)<10e-3:
for i=1:length(vec_delx_out)
    if norm(x_vec_out(i,:)) < 10e-3
        break
    end
end

norm_x = [];
for j=1:length(t_out)
    norm_x(j) = norm(x_vec_out(j,:));
end

fig2 = figure ("Name","2D Figure of The Orbit and Target Trajectory Over Time",'Position',[300 300 900 500]);
hold all

plot(t_out, x_vec_out(:,1), "LineWidth", 2, "Color", colors(1,:))
plot(t_out, x_vec_out(:,2), "LineWidth", 2, "Color", colors(3,:))
plot(t_out, x_vec_out(:,3), "LineWidth", 2, "Color", colors(5,:))
plot(t_out, norm_x, "LineWidth", 2, "Color", colors(7,:))
plot(t_out, x_ref_vec_out(:,1),"-", "LineWidth", 2, "Color", "r")
plot(t_out, x_ref_vec_out(:,2),"--", "LineWidth", 2, "Color", "k")
plot(t_out, x_ref_vec_out(:,3),"--", "LineWidth", 2, "Color", "b")
plot(ones(length(linspace(0,0.5,50)))*t_out(i), linspace(0,0.5,50), "LineWidth", 2, "Color", colors(9,:))

xlabel('Time [sec]')
ylabel('Position [km]')
grid on
grid minor
title("2D Figure of The Orbit and Target Trajectory Over Time")
subtitle("Almog Dobrescu 214254252")
legend({'x', '$x_{ref}$', 'y', '$y_{ref}$', 'z', '$z_{ref}$', 'distance', '$|x|<10$'},'FontSize',11 ,'Location','northeast',Interpreter='latex')
% exportgraphics(fig2, 'graph2.png','Resolution',600);
%%
fig3 = figure ("Name","Thrust Acceleration Components and Total Thrust Acceleration Over Time",'Position',[500 300 900 500]);
colors = cool(4);
fs = zeros(length(t_out),3);
for j = 1:length(t_out)
    t = t_out(j);
    delx_vec_state = state(1:6, :);
    f = K*delx_vec_state;
    norm_f = norm(f);
    if norm_f > f_max
        f = f/norm_f*f_max*0.9;
    end
    fs(j,:) = f;
end

hold all

plot(t_out, fs(:,1), "LineWidth", 2, "Color", colors(1,:))
plot(t_out, fs(:,2), "LineWidth", 2, "Color", colors(2,:))
plot(t_out, fs(:,3), "LineWidth", 2, "Color", colors(3,:))
plot(t_out, sqrt(fs(:,1).^2+fs(:,2).^2+fs(:,3).^2), "LineWidth", 2, "Color", colors(4,:))
plot(t_out, ones(length(t_out),1)*f_max,"--", "LineWidth", 2, "Color", "k")
plot(t_out, -ones(length(t_out),1)*f_max,"--", "LineWidth", 2, "Color", "k")

xlabel('Time [sec]','FontSize', 16, 'Interpreter','latex')
ylabel('Acceleration $\left[\frac{km}{sec^2}\right]$','FontSize', 16, 'Interpreter','latex')
grid on
grid minor
title("Thrust Acceleration Components and Total Thrust Acceleration Over Time")
subtitle("Almog Dobrescu 214254252")
legend({'f_x', 'f_y', 'f_z', 'f_t_o_t', 'f_l_i_m_i_t'},'FontSize',11 ,'Location','northeast')
% exportgraphics(fig3, 'graph3.png','Resolution',600);
%%
fig4 = figure ("Name","Total delta v Over Time",'Position',[700 300 900 500]);

delta_v = [];
for j = 2:length(t_out)
    delta_v(j) = trapz(t_out(1:j), sqrt(fs(1:j,1).^2+fs(1:j,2).^2+fs(1:j,3).^2));
end
hold all

plot([t_out], [delta_v], "LineWidth", 2, "Color", "k")

xlabel('Time [sec]','FontSize', 16, 'Interpreter','latex')
ylabel('$\Delta v\left[\frac{km}{sec}\right]$','FontSize', 16, 'Interpreter','latex')
grid on
grid minor
title("Total $\Delta v$ Over Time","Interpreter","latex")
subtitle("Almog Dobrescu 214254252")
% legend({''},'FontSize',11 ,'Location','northeast')
% exportgraphics(fig4, 'graph4.png','Resolution',600);
%%
colors = cool(2);

for j=1:length(t_out)
    norm_delx(j) = norm(vec_delx_out(j,:));
end

fig5 = figure ("Name","The Miss Distance From the Target Trajectory Over Time",'Position',[400 100 900 500]);
hold all

plot(t_out, norm_delx, "LineWidth", 2, "Color", colors(1,:))
plot(ones(length(linspace(0,0.5,50)))*t_out(i), linspace(0,0.5,50), "LineWidth", 2, "Color", colors(2,:))

xlabel("$Time [sec]$",Interpreter="latex")
ylabel("$\delta x [km]$",Interpreter="latex")
grid on
grid minor
title("The Miss Distance From the Target Trajectory Over Time")
subtitle("Almog Dobrescu 214254252")
legend({'$\left|\delta x\right|$', '$\left|x\right|<10$'},'FontSize',11 ,'Location','northeast',Interpreter='latex')
% exportgraphics(fig5, 'graph5.png','Resolution',600);
