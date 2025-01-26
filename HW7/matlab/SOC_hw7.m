%% A
clc; clear; close all;

global_variabels

state = [x0_vec;x0_vec;x0_vec;x0_vec]; % [x3;v3;x4;v4;x6;v6;x7;v7]
time_interval = [0:t_f];

 % This is where we integrate the equations of motion.
[t_out, state_out] = ode45(@Satellite, time_interval, state, odeset('RelTol',5e-14,'AbsTol',5e-14));

x3_vec_out = [state_out(:,1),state_out(:,3),state_out(:,5)];
v3_vec_out = [state_out(:,2),state_out(:,4),state_out(:,6)];
x4_vec_out = [state_out(:,7),state_out(:,9),state_out(:,11)];
v4_vec_out = [state_out(:,8),state_out(:,10),state_out(:,12)];
x6_vec_out = [state_out(:,13),state_out(:,15),state_out(:,17)];
v6_vec_out = [state_out(:,14),state_out(:,16),state_out(:,18)];
x7_vec_out = [state_out(:,19),state_out(:,21),state_out(:,23)];
v7_vec_out = [state_out(:,20),state_out(:,22),state_out(:,24)];

miss_distance = norm([x7_vec_out(end,1),x7_vec_out(end,2),x7_vec_out(end,3)]);
miss_velocity = norm([v7_vec_out(end,1),v7_vec_out(end,2),v7_vec_out(end,3)]);

%%
fig1 = figure ("Name","3D Figure of The Orbit Trajectory",'Position',[100 300 900 500]);
colors = cool(4)*0.9;
hold all
plot3(x7_vec_out(:,1),x7_vec_out(:,2),x7_vec_out(:,3), "LineWidth", 2, "Color", colors(1,:))
plot3(x6_vec_out(:,1),x6_vec_out(:,2),x6_vec_out(:,3), "LineWidth", 2, "Color", colors(2,:))
plot3(x4_vec_out(:,1),x4_vec_out(:,2),x4_vec_out(:,3), "LineWidth", 2, "Color", colors(3,:))
plot3(x3_vec_out(:,1),x3_vec_out(:,2),x3_vec_out(:,3), "LineWidth", 2, "Color", colors(4,:))

plot3(x7_vec_out(1,1),x7_vec_out(1,2),x7_vec_out(1,3), "o", "LineWidth", 4, "Color", "k")
plot3(x7_vec_out(end,1),x7_vec_out(end,2),x7_vec_out(end,3),'pentagram', "LineWidth", 2, "Color", "r")

xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
grid on
grid minor
title(sprintf("3D Figure of The Orbit Trajectory | miss distance = %g[m] | miss velocity = %g[cm/sec]", miss_distance*10^3, miss_velocity*10^5))
subtitle("Almog Dobrescu 214254252")
legend({'Linear Time-Varying','LQR', 'Pole Placement', 'Cross Product', 'init position', 'end position'},'FontSize',11 ,'Location','northwest')
% exportgraphics(fig1, 'graph1.png','Resolution',300);
%%
linear_colors = copper(4)*0.9;
lqr_colors = cool(4)*0.9;
pp_colors = summer(4)*0.9;
cp_colors = winter(4)*0.9;

norm_x3 = [];
norm_x4 = [];
norm_x6 = [];
norm_x7 = [];
for j=1:length(t_out)
    norm_x3(j) = norm(x3_vec_out(j,:));
    norm_x4(j) = norm(x4_vec_out(j,:));
    norm_x6(j) = norm(x6_vec_out(j,:));
    norm_x7(j) = norm(x7_vec_out(j,:));
end

fig2 = figure ("Name","2D Figure of The Orbit",'Position',[300 300 900 500]);
hold all
% 7
plot(t_out, x7_vec_out(:,1),"-", "LineWidth", 3, "Color", linear_colors(1,:))
plot(t_out, x7_vec_out(:,2),"-", "LineWidth", 3, "Color", linear_colors(2,:))
plot(t_out, x7_vec_out(:,3),"-", "LineWidth", 3, "Color", linear_colors(3,:))
plot(t_out, norm_x7,"-", "LineWidth", 3, "Color", linear_colors(4,:))
% 6
plot(t_out, x6_vec_out(:,1),"-", "LineWidth", 1, "Color", lqr_colors(1,:))
plot(t_out, x6_vec_out(:,2),"-", "LineWidth", 1, "Color", lqr_colors(2,:))
plot(t_out, x6_vec_out(:,3),"-", "LineWidth", 1, "Color", lqr_colors(3,:))
plot(t_out, norm_x6,"-", "LineWidth", 1, "Color", lqr_colors(4,:))
% 4
plot(t_out, x4_vec_out(:,1),"--", "LineWidth", 1.5, "Color", pp_colors(1,:))
plot(t_out, x4_vec_out(:,2),"--", "LineWidth", 1.5, "Color", pp_colors(2,:))
plot(t_out, x4_vec_out(:,3),"--", "LineWidth", 1.5, "Color", pp_colors(3,:))
plot(t_out, norm_x4,"--", "LineWidth", 1.5, "Color", pp_colors(4,:))
% 3
plot(t_out, x3_vec_out(:,1),"-.", "LineWidth", 1.5, "Color", cp_colors(1,:))
plot(t_out, x3_vec_out(:,2),"-.", "LineWidth", 1.5, "Color", cp_colors(2,:))
plot(t_out, x3_vec_out(:,3),"-.", "LineWidth", 1.5, "Color", cp_colors(3,:))
plot(t_out, norm_x3,"-.", "LineWidth", 1.5, "Color", cp_colors(4,:))

xlabel('Time [sec]')
ylabel('Position [km]')
grid on
grid minor
title("2D Figure of The Orbit Trajectory Over Time")
subtitle("Almog Dobrescu 214254252")
legend({'Linear Time-Varying-x', 'Linear Time-Varying-y', 'Linear Time-Varying-z', 'Linear Time-Varying-norm', 'LQR-x', 'LQR-y', 'LQR-z', 'LQR-norm', 'Pole Placement-x', 'Pole Placement-y', 'Pole Placement-z', 'Pole Placement-norm', 'Cross Product-x', 'Cross Product-y', 'Cross Product-z', 'Cross Product-norm'},'FontSize',11 ,'Location','northeast')
% exportgraphics(fig2, 'graph2.png','Resolution',300);
%%
fig3 = figure ("Name","Thrust Acceleration Components and Total Thrust Acceleration Over Time",'Position',[500 300 900 500]);
linear_colors = copper(4)*0.9;
lqr_colors = cool(4)*0.9;
pp_colors = summer(4)*0.9;
cp_colors = winter(4)*0.9;

fs3 = zeros(length(t_out),3);
fs4 = zeros(length(t_out),3);
fs6 = zeros(length(t_out),3);
fs7 = zeros(length(t_out),3);

for j = 1:length(t_out)
    t = t_out(j);
    x3_vec = state_out(j, 1:6);
    pos3_vec = [x3_vec(1);x3_vec(3);x3_vec(5)];
    v3_vec = [x3_vec(2);x3_vec(4);x3_vec(6)];
    x4_vec = state_out(j, 7:12);
    x6_vec = state_out(j, 13:18);
    
    %3
    C_star = calc_C_star(t_f, t, n);
    [phi_11, phi_12, ~, ~] = calc_phis(t_f, t, n);
    vg3 = C_star*pos3_vec - v3_vec;
    vg3_normaliz = vg3./norm(vg3);
    p3 = -C_star*vg3;
    q3 = (f_max^2 - norm(p3)^2 + (dot(p3,vg3_normaliz))^2)^0.5;
    f3 = p3 + (q3 - dot(p3,vg3_normaliz))*vg3_normaliz;
    if norm(phi_11*pos3_vec + phi_12*v3_vec) < 1e-3
        f3 = zeros(3,1);
    end
    if ~isreal(q3) || ~isfinite(q3)
        f3 = zeros(3,1);
    end
    fs3(j,:) = f3;

    %4
    f4 = -K4*x4_vec.';
    if norm(f4) > f_max
        f4 = f4/norm(f4)*f_max;
    end
    fs4(j,:) = f4;

    %6
    f6 = -K6*x6_vec.';
    norm_f6 = norm(f6);
    if norm_f6 > f_max
        f6 = f6/norm_f6*f_max;
    end
    fs6(j,:) = f6;

    %7
    fun=@(tau) expm(F*(t_f-tau))*(G*G')*expm(-F'*tau);
    M=-integral(fun,0,t_f,'ArrayValued',true);
    nu=-(Psi*M*expm(F'*t_f)*Psi')^-1*Psi*expm(F*t_f)*x0_vec;
    lambda = expm(-F'*(t-t_f))*Psi'*nu;
    f7 = -G'*lambda;
    
    if ~isfinite(f7)
        f7 = zeros(3,1);
    end
    norm_f7 = norm(f7);
    if norm_f7 > f_max
        f7 = f7/norm_f7*f_max;
    end
    fs7(j,:) = f7;

end

hold all

%7
plot(t_out, fs7(:,1), "-", "LineWidth", 3, "Color", linear_colors(1,:))
plot(t_out, fs7(:,2), "-", "LineWidth", 3, "Color", linear_colors(2,:))
plot(t_out, fs7(:,3), "-", "LineWidth", 3, "Color", linear_colors(3,:))
plot(t_out, sqrt(fs7(:,1).^2+fs7(:,2).^2+fs7(:,3).^2), "-", "LineWidth", 3, "Color", linear_colors(4,:))
%6
plot(t_out, fs6(:,1), "-", "LineWidth", 1.5, "Color", lqr_colors(1,:))
plot(t_out, fs6(:,2), "-", "LineWidth", 1.5, "Color", lqr_colors(2,:))
plot(t_out, fs6(:,3), "-", "LineWidth", 1.5, "Color", lqr_colors(3,:))
plot(t_out, sqrt(fs6(:,1).^2+fs6(:,2).^2+fs6(:,3).^2), "-", "LineWidth", 1.5, "Color", lqr_colors(4,:))
%4
plot(t_out, fs4(:,1), "--", "LineWidth", 1.5, "Color", pp_colors(1,:))
plot(t_out, fs4(:,2), "--", "LineWidth", 1.5, "Color", pp_colors(2,:))
plot(t_out, fs4(:,3), "--", "LineWidth", 1.5, "Color", pp_colors(3,:))
plot(t_out, sqrt(fs4(:,1).^2+fs4(:,2).^2+fs4(:,3).^2), "--", "LineWidth", 1.5, "Color", pp_colors(4,:))
%3
plot(t_out, fs3(:,1), "-.", "LineWidth", 1.5, "Color", cp_colors(1,:))
plot(t_out, fs3(:,2), "-.", "LineWidth", 1.5, "Color", cp_colors(2,:))
plot(t_out, fs3(:,3), "-.", "LineWidth", 1.5, "Color", cp_colors(3,:))
plot(t_out, sqrt(fs3(:,1).^2+fs3(:,2).^2+fs3(:,3).^2), "-.", "LineWidth", 1.5, "Color", cp_colors(4,:))

plot(t_out, ones(length(t_out),1)*f_max,":", "LineWidth", 2, "Color", "k")
plot(t_out, -ones(length(t_out),1)*f_max,":", "LineWidth", 2, "Color", "k")

xlabel('Time [sec]','FontSize', 16, 'Interpreter','latex')
ylabel('Acceleration $\left[\frac{km}{sec^2}\right]$','FontSize', 16, 'Interpreter','latex')
grid on
grid minor
title("Thrust Acceleration Components and Total Thrust Acceleration Over Time")
subtitle("Almog Dobrescu 214254252")
legend({'Linear Time-Varying-fx', 'Linear Time-Varying-fy', 'Linear Time-Varying-fz', 'Linear Time-Varying-f_norm','LQR-fx', 'LQR-fy', 'LQR-fz', 'LQR-f_norm', 'Pole Placement-fx', 'Pole Placement-fy', 'Pole Placement-fz', 'Pole Placement-f_norm', 'Cross Product-fx', 'Cross Product-fy', 'Cross Product-fz', 'Cross Product-f_norm', 'f_max'},'FontSize',11 ,'Location','northeast',Interpreter='none')
% exportgraphics(fig3, 'graph3.png','Resolution',300);
%%
fig4 = figure ("Name","Thrust Acceleration Components and Total Thrust Acceleration Over Time",'Position',[500 300 900 500]);
linear_colors = copper(4)*0.9;

hold all

%7
plot(t_out, fs7(:,1), "-", "LineWidth", 1.5, "Color", linear_colors(1,:))
plot(t_out, fs7(:,2), "-", "LineWidth", 1.5, "Color", linear_colors(2,:))
plot(t_out, fs7(:,3), "-", "LineWidth", 1.5, "Color", linear_colors(3,:))
plot(t_out, sqrt(fs7(:,1).^2+fs7(:,2).^2+fs7(:,3).^2), "-", "LineWidth", 1.5, "Color", linear_colors(4,:))

% plot(t_out, ones(length(t_out),1)*f_max,":", "LineWidth", 2, "Color", "k")
% plot(t_out, -ones(length(t_out),1)*f_max,":", "LineWidth", 2, "Color", "k")

xlabel('Time [sec]','FontSize', 16, 'Interpreter','latex')
ylabel('Acceleration $\left[\frac{km}{sec^2}\right]$','FontSize', 16, 'Interpreter','latex')
grid on
grid minor
title("Thrust Acceleration Components and Total Thrust Acceleration Over Time")
subtitle("Almog Dobrescu 214254252")
legend({'Linear Time-Varying-fx', 'Linear Time-Varying-fy', 'Linear Time-Varying-fz', 'Linear Time-Varying-f_norm', 'f_max'},'FontSize',11 ,'Location','northeast',Interpreter='none')
% exportgraphics(fig4, 'graph4.png','Resolution',300);

%%
fig5 = figure ("Name","Total delta v Over Time",'Position',[700 300 900 500]);
linear_colors = copper(4)*0.9;
lqr_colors = cool(4)*0.9;
pp_colors = summer(4)*0.9;
cp_colors = winter(4)*0.9;

delta_v3 = [];
delta_v4 = [];
delta_v6 = [];
delta_v7 = [];
for j = 2:length(t_out)
    delta_v3(j) = trapz(t_out(1:j), sqrt(fs3(1:j,1).^2+fs3(1:j,2).^2+fs3(1:j,3).^2));
    delta_v4(j) = trapz(t_out(1:j), sqrt(fs4(1:j,1).^2+fs4(1:j,2).^2+fs4(1:j,3).^2));
    delta_v6(j) = trapz(t_out(1:j), sqrt(fs6(1:j,1).^2+fs6(1:j,2).^2+fs6(1:j,3).^2));
    delta_v7(j) = trapz(t_out(1:j), sqrt(fs7(1:j,1).^2+fs7(1:j,2).^2+fs7(1:j,3).^2));
end
hold all

plot([t_out], [delta_v7], "LineWidth", 1.5, "Color", linear_colors(1,:))
plot([t_out], [delta_v6], "LineWidth", 1.5, "Color", lqr_colors(1,:))
plot([t_out], [delta_v4], "LineWidth", 1.5, "Color", pp_colors(1,:))
plot([t_out], [delta_v3], "LineWidth", 1.5, "Color", cp_colors(1,:))


xlabel('Time [sec]','FontSize', 16, 'Interpreter','latex')
ylabel('$\Delta v\left[\frac{km}{sec}\right]$','FontSize', 16, 'Interpreter','latex')
grid on
grid minor
title("Total $\Delta v$ Over Time","Interpreter","latex")
subtitle("Almog Dobrescu 214254252")
legend({'Linear Time-Varying-$\Delta v$', 'LQR-$\Delta v$', 'Pole Placement-$\Delta v$', 'Cross Product-$\Delta v$'},'FontSize',11 ,'Location','northeast', "Interpreter","latex")
% exportgraphics(fig5, 'graph5.png','Resolution',300);
