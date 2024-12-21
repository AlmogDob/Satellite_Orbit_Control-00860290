function d_state_dt = satellite_caseB(t, state)
x_vec = state(1:6, 1);
pos_vec = [state(1);state(3);state(5)];
v_vec = [state(2);state(4);state(6)];

global_variabels
t
tau = t_f - t;
vg = double(subs(C_star)*pos_vec - v_vec);

if norm(vg) < 1e-8
    vg = zeros(3,1);
end

p = double(-subs(C_star)*vg);
q = (accel_limit^2 - norm(p)^2 + (dot(p,vg/norm(vg)))^2)^0.5;

f = p + (q - dot(p,vg/norm(vg)))*vg/norm(vg);
if norm(f)>accel_limit
    f = f/norm(f)*accel_limit*0.9;
end


d_x_dt  = F*x_vec + G*f;

d_state_dt = [d_x_dt];

