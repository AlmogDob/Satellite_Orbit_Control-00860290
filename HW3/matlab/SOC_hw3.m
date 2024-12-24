%% A
clc; clear;

global_variabels

pos0_vec = [0;-1;1]; %[km]
v0_vec = [0;0;-0.74267*n]; %[km/sec]
posf_vec = [0;0;0];
vf_vec = [0;0;0];
v0_req_vec = (posf_vec-pos0_vec)./t_f;

delta_v1 = v0_req_vec - v0_vec;
delta_vf = vf_vec - v0_req_vec;

max_f = sqrt((2*n*0.5e-3)^2+(n^2)^2);

state = [pos0_vec(1);v0_req_vec(1);pos0_vec(2);v0_req_vec(2);pos0_vec(3);v0_req_vec(3)];
time_interval = [0, t_f];

 % This is where we integrate the equations of motion.
[t_out, state_out] = ode45(@satellite_caseA, time_interval, state, odeset('RelTol',1e-10,'AbsTol',1e-10));

miss_distance = norm([state_out(end,1),state_out(end,3),state_out(end,5)]);
miss_velocity = norm([state_out(end,2),state_out(end,4),state_out(end,6)]);

fig1 = figure ("Name","3D Figure of The Orbit Trajectory",'Position',[100 300 900 500]);
hold all

plot3(state_out(:,1),state_out(:,3),state_out(:,5), "LineWidth", 2, "Color", "#D95319")
plot3(state_out(1,1),state_out(1,3),state_out(1,5), "o", "LineWidth", 4, "Color", "k")
plot3(state_out(end,1),state_out(end,3),state_out(end,5),'hexagram', "LineWidth", 4, "Color", "m")

xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
grid on
grid minor
title(sprintf("3D Figure of The Orbit Trajectory | miss distance = %g[m]", miss_distance*10^3))
subtitle("Almog Dobrescu 214254252")
legend({'the orbit', 'init position', 'end position'},'FontSize',11 ,'Location','northwest')
% exportgraphics(fig1, 'graph1.png','Resolution',600);

fig2 = figure ("Name","2D Figure of The Orbit Trajectory Over Time",'Position',[300 300 900 500]);
hold all

plot(t_out, state_out(:,1), "LineWidth", 2)
plot(t_out, state_out(:,3), "LineWidth", 2)
plot(t_out, state_out(:,5), "LineWidth", 2)

xlabel('Time [sec]')
ylabel('Position [km]')
grid on
grid minor
title("2D Figure of The Orbit Trajectory Over Time")
subtitle("Almog Dobrescu 214254252")
legend({'x', 'y', 'z'},'FontSize',11 ,'Location','northeast')
% exportgraphics(fig2, 'graph2.png','Resolution',600);

fig3 = figure ("Name","Thrust Acceleration Components and Total Thrust Acceleration Over Time",'Position',[500 300 900 500]);
colors = cool(4);
f = [-2*n*0.5e-3*ones(length(t_out),1), zeros(length(t_out),1), n^2.*state_out(:,5)];
hold all

plot(t_out, f(:,1), "LineWidth", 2, "Color", colors(1,:))
plot(t_out, f(:,2), "LineWidth", 2, "Color", colors(2,:))
plot(t_out, f(:,3), "LineWidth", 2, "Color", colors(3,:))
plot(t_out, sqrt(f(:,1).^2+f(:,2).^2+f(:,3).^2), "LineWidth", 2, "Color", colors(4,:))
plot(t_out, ones(length(t_out),1)*accel_limit,"--", "LineWidth", 2, "Color", "k")
plot(t_out, -ones(length(t_out),1)*accel_limit,"--", "LineWidth", 2, "Color", "k")

xlabel('Time [sec]','FontSize', 16, 'Interpreter','latex')
ylabel('Acceleration $\left[\frac{km}{sec^2}\right]$','FontSize', 16, 'Interpreter','latex')
grid on
grid minor
title("Thrust Acceleration Components and Total Thrust Acceleration Over Time")
subtitle("Almog Dobrescu 214254252")
legend({'f_x', 'f_y', 'f_z', 'f_t_o_t', 'f_l_i_m_i_t'},'FontSize',11 ,'Location','northeast')
% exportgraphics(fig3, 'graph3.png','Resolution',600);

fig4 = figure ("Name","Total delta v Over Time",'Position',[700 300 900 500]);
f = [-2*n*0.5e-3*ones(length(t_out),1), zeros(length(t_out),1), n^2.*state_out(:,5)];
delta_v = [];
for i = 2:length(t_out)
    delta_v(i) = trapz(t_out(1:i), sqrt(f(1:i,1).^2+f(1:i,2).^2+f(1:i,3).^2));
end
hold all

plot([0;0;t_out;t_out(end);t_out(end)], [0, norm(delta_v1), norm(delta_v1)+delta_v, norm(delta_v1)+delta_v(end), norm(delta_v1)+delta_v(end)+norm(delta_vf)], "LineWidth", 2, "Color", "k")

xlabel('Time [sec]','FontSize', 16, 'Interpreter','latex')
ylabel('$\Delta v\left[\frac{km}{sec}\right]$','FontSize', 16, 'Interpreter','latex')
grid on
grid minor
title("Total $\Delta v$ Over Time","Interpreter","latex")
subtitle("Almog Dobrescu 214254252")
% legend({''},'FontSize',11 ,'Location','northeast')
% exportgraphics(fig4, 'graph4.png','Resolution',600);


%% B
clc; clear;

global_variabels

state = [pos0_vec(1);v0_vec(1);pos0_vec(2);v0_vec(2);pos0_vec(3);v0_vec(3)];
time_interval = [0, t_f-1e-5];
tol = 1e-10;
 % This is where we integrate the equations of motion.
[t_out, state_out] = ode45(@satellite_caseB, time_interval, state, odeset('RelTol',tol,'AbsTol',tol));

miss_distance = norm([state_out(end,1),state_out(end,3),state_out(end,5)]);
miss_velocity = norm([state_out(end,2),state_out(end,4),state_out(end,6)]);

fig5 = figure ("Name","3D Figure of The Orbit Trajectory",'Position',[100 150 900 500]);
hold all

plot3(state_out(:,1),state_out(:,3),state_out(:,5), "LineWidth", 2, "Color", "#D95319")
plot3(state_out(1,1),state_out(1,3),state_out(1,5), "o", "LineWidth", 4, "Color", "k")
plot3(state_out(end,1),state_out(end,3),state_out(end,5),'hexagram', "LineWidth", 4, "Color", "m")

xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
grid on
grid minor
title(sprintf("3D Figure of The Orbit Trajectory | miss distance = %g[m] | miss velocity = %g[m/sec]", miss_distance*10^3, miss_velocity*10^3))
subtitle("Almog Dobrescu 214254252")
legend({'the orbit', 'init position', 'end position'},'FontSize',11 ,'Location','northwest')
% exportgraphics(fig5, 'graph5.png','Resolution',600);
%%
fig6 = figure ("Name","2D Figure of The Orbit Trajectory Over Time",'Position',[300 150 900 500]);
hold all

plot(t_out, state_out(:,1), "LineWidth", 2)
plot(t_out, state_out(:,3), "LineWidth", 2)
plot(t_out, state_out(:,5), "LineWidth", 2)

xlabel('Time [sec]')
ylabel('Position [km]')
grid on
grid minor
title("2D Figure of The Orbit Trajectory Over Time")
subtitle("Almog Dobrescu 214254252")
legend({'x', 'y', 'z'},'FontSize',11 ,'Location','northeast')
% exportgraphics(fig6, 'graph6.png','Resolution',600);
%%
fig7 = figure ("Name","Thrust Acceleration Components and Total Thrust Acceleration Over Time",'Position',[500 150 900 500]);
colors = cool(4);
fs = zeros(length(t_out),3);
for i = 1:length(t_out)
    t = t_out(i);
    pos_vec = [state_out(i,1);state_out(i,3);state_out(i,5)];
    v_vec   = [state_out(i,2);state_out(i,4);state_out(i,6)];
    C_star = calc_C_star(t_f, t, n);
    [phi_11, phi_12, ~, ~] = calc_phis(t_f, t, n);
    
    vg = C_star*pos_vec - v_vec;
    
    vg_normaliz = vg./norm(vg);
    
    p = -C_star*vg;
    q = (accel_limit^2 - norm(p)^2 + (dot(p,vg_normaliz))^2)^0.5;
    
    f = p + (q - dot(p,vg_normaliz))*vg_normaliz;
    % if norm(vg) < 1e-3
    %     f = zeros(3,1);
    % end
    if norm(phi_11*pos_vec + phi_12*v_vec) < 1e-3
        f = zeros(3,1);
    end
    if ~isreal(q)
        f = zeros(3,1);
    end
    fs(i,:) = f;
end

hold all

plot(t_out, fs(:,1), "LineWidth", 2, "Color", colors(1,:))
plot(t_out, fs(:,2), "LineWidth", 2, "Color", colors(2,:))
plot(t_out, fs(:,3), "LineWidth", 2, "Color", colors(3,:))
plot(t_out, sqrt(fs(:,1).^2+fs(:,2).^2+fs(:,3).^2), "LineWidth", 2, "Color", colors(4,:))
plot(t_out, ones(length(t_out),1)*accel_limit,"--", "LineWidth", 2, "Color", "k")
plot(t_out, -ones(length(t_out),1)*accel_limit,"--", "LineWidth", 2, "Color", "k")

xlabel('Time [sec]','FontSize', 16, 'Interpreter','latex')
ylabel('Acceleration $\left[\frac{km}{sec^2}\right]$','FontSize', 16, 'Interpreter','latex')
grid on
grid minor
title("Thrust Acceleration Components and Total Thrust Acceleration Over Time")
subtitle("Almog Dobrescu 214254252")
legend({'f_x', 'f_y', 'f_z', 'f_t_o_t', 'f_l_i_m_i_t'},'FontSize',11 ,'Location','northeast')
% exportgraphics(fig7, 'graph7.png','Resolution',600);
%%
fig8 = figure ("Name","Total delta v Over Time",'Position',[700 150 900 500]);

delta_v = [];
for i = 2:length(t_out)
    delta_v(i) = trapz(t_out(1:i), sqrt(fs(1:i,1).^2+fs(1:i,2).^2+fs(1:i,3).^2));
end
hold all

plot([t_out;t_out(end);t_out(end)], [delta_v, delta_v(end), delta_v(end)+norm(miss_velocity)], "LineWidth", 2, "Color", "k")

xlabel('Time [sec]','FontSize', 16, 'Interpreter','latex')
ylabel('$\Delta v\left[\frac{km}{sec}\right]$','FontSize', 16, 'Interpreter','latex')
grid on
grid minor
title("Total $\Delta v$ Over Time","Interpreter","latex")
subtitle("Almog Dobrescu 214254252")
% legend({''},'FontSize',11 ,'Location','northeast')
% exportgraphics(fig8, 'graph8.png','Resolution',600);


