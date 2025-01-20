function d_state_dt = Satellite(t, state)
x3_vec = state(1:6, :);
pos3_vec = [x3_vec(1);x3_vec(3);x3_vec(5)];
v3_vec = [x3_vec(2);x3_vec(4);x3_vec(6)];

x4_vec = state(7:12, :);
x6_vec = state(13:18, :);

global_variabels

%% 3
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
d_x3_dt  = F*x3_vec + G*f3;

%% 4
f4 = -K4*x4_vec;
if norm(f4) > f_max
    f4 = f4/norm(f4)*f_max;
end

d_x4_dt = F*x4_vec + G*f4;

%% 6
f6 = -K6*x6_vec;
norm_f6 = norm(f6);
if norm_f6 > f_max
    f6 = f6/norm_f6*f_max;
end
d_x6_dt = F*x6_vec + G*f6;

d_state_dt = [d_x3_dt;d_x4_dt;d_x6_dt];

