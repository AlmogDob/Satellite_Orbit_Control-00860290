function d_state_dt = Satellite(t, state)
delx_vec = state(1:6, :);
xr_vec = state(7:12, :);

global_variabels

f = K*delx_vec;
norm_f = norm(f);
if norm_f > f_max
    f = f/norm_f*f_max*0.9;
end

Gf      = -(F-Fr)*xr_vec - G*f;
d_delx_dt  = F*delx_vec + (F-Fr)*xr_vec + Gf;
d_xr_dt = Fr*xr_vec;

d_state_dt = [d_delx_dt; d_xr_dt];

