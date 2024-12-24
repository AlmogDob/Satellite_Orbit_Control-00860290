function [C_star] = calc_C_star(t_f, t, n)
[phi_11, phi_12, ~, ~] = calc_phis(t_f, t, n);
C_star = -inv(phi_12)*phi_11;
end