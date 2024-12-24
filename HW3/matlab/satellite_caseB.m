function d_state_dt = satellite_caseB(t, state)
x_vec = state(1:6, 1);
pos_vec = [state(1);state(3);state(5)];
v_vec = [state(2);state(4);state(6)];

global_variabels

t
C_star = calc_C_star(t_f, t, n);
[phi_11, phi_12, ~, ~] = calc_phis(t_f, t, n);

vg = C_star*pos_vec - v_vec;

vg_normaliz = vg./norm(vg);

p = -C_star*vg;
q = (accel_limit^2 - norm(p)^2 + (dot(p,vg_normaliz))^2)^0.5;

f = p + (q - dot(p,vg_normaliz))*vg_normaliz;
% if norm(vg) < 1e-4
%     f = zeros(3,1);
% end
if norm(phi_11*pos_vec + phi_12*v_vec) < 1e-3
    f = zeros(3,1);
end
if ~isreal(q)
    f = zeros(3,1);
end
d_x_dt  = F*x_vec + G*f;

d_state_dt = [d_x_dt];
end
