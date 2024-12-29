function d_state_dt = SatelliteC(t, state)
x_vec = state(1:6, 1);
global_variabels

if t > 500
    x_vec(5, 1) = 0;
end

f = -K*x_vec;
if norm(f) > f_max
    f = f/norm(f)*f_max;
end

d_x_dt  = F*x_vec + G*f;

d_state_dt = [d_x_dt];

