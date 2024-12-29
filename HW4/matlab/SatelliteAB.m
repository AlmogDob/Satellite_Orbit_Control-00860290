function d_state_dt = SatelliteAB(t, state)

global_variabels

x_vec = state(1:6, :);
t;

f = -K*x_vec;
if norm(f) > f_max
    f = f/norm(f)*f_max;
end

d_x_dt = F*x_vec + G*f;

d_state_dt = [d_x_dt];

