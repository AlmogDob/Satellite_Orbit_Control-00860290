function d_state_dt = Satellite(t, state)
x_vec = state(1:6, :);

global_variabels

f = [-2*n*0.5e-3;
          0     ;
     n^2*x_vec(5)];

d_x_dt  = F*x_vec + G*f;

d_state_dt = [d_x_dt];

