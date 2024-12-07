function d_state_dt = Satellite(t, state)
global_variabels;

r_vector = state(1:3);
v_vector = state(4:6);

f_drag = calculate_drag_force(r_vector, v_vector);
f_j2 = calculate_J2_force(r_vector);
f_s = calculate_sum_pressure_force(r_vector);

r = norm(r_vector);
r_dot_vector = v_vector;
% v_dot_vector = -myu*r_vector/r^3 + f_drag - f_j2;% + f_s;
% v_dot_vector = -myu*r_vector/r^3 + f_drag + f_s;
v_dot_vector = -myu*r_vector/r^3 + f_drag - f_j2 + f_s;


d_state_dt = [r_dot_vector;v_dot_vector];
