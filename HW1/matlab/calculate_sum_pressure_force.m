function [f] = calculate_sum_pressure_force(r_vec)
global_variabels
in_the_shadow = 0;

N = datenum(['1-Jan']) - datenum(['1-Jan']) + 1;
RA_deg = 0.98563*(N-80);
delta_deg = asind(0.39795 * cosd(0.98563*(N - 173)));

r_s_vec = [cosd(delta_deg)*cosd(RA_deg);
           cosd(delta_deg)*sind(RA_deg);
           sind(delta_deg)];

if (dot(r_s_vec, r_vec) < 0)
    in_the_shadow = 1;
end
if (norm(r_vec) * sin(acos(dot(r_vec, r_s_vec))) < R_earth)
    in_the_shadow = 1; 
end

if (in_the_shadow)
    f = 0;
else
    f = -(1+epsilon) * G_s/c * A_solar/m * r_s_vec;
end