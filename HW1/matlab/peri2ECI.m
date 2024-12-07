function [r_ECI] = peri2ECI(vector_per, i, small_omega, big_omega)

first_matrix_rotation = [cos(big_omega) sin(big_omega) 0;-sin(big_omega) cos(big_omega) 0;0 0 1];
second_matrix_rotation = [1 0 0;0 cos(i) sin(i);0 -sin(i) cos(i)];
third_matrix_rotation = [cos(small_omega) sin(small_omega) 0;-sin(small_omega) cos(small_omega) 0;0 0 1];
ECI2peri = third_matrix_rotation * second_matrix_rotation * first_matrix_rotation;
r_ECI = transpose(ECI2peri) * vector_per;
end