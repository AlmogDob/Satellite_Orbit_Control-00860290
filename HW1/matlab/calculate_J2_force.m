function [f] = calculate_J2_force(r_vec)
global_variabels

x = r_vec(1);
y = r_vec(2);
z = r_vec(3);
r = norm(r_vec);

f = zeros(3,1);

f(1) = -3/2 * myu*x/r^3 * J_2 * (R_earth/r)^2 * (5*z^2/r^2 - 1);
f(2) = -3/2 * myu*y/r^3 * J_2 * (R_earth/r)^2 * (5*z^2/r^2 - 1);
f(3) = -3/2 * myu*z/r^3 * J_2 * (R_earth/r)^2 * (5*z^2/r^2 - 3);
end