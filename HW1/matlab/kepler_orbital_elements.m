function [a,e,i,small_omega, big_omega,f] = kepler_orbital_elements(r_vector,v_vector, mu)
%kepler_orbital_elements retrurns the kepler orbiral elements using the
%position and the velocity vectors
h_vector = cross(r_vector, v_vector);
h = norm(h_vector);
r = norm(r_vector);
v = norm(v_vector);
e_vector = cross(v_vector, h_vector)/mu - r_vector/r;
e = norm(e_vector);
p = h^2/mu;
a = p/(1-e^2);
n_vector = cross([0 0 1], h_vector);
n = norm(n_vector);
i = acos(h_vector(3)/h);
big_omega = atan2(n_vector(2)/n, n_vector(1)/n);
small_omega = atan2(sign(e_vector(3))*sqrt(1-(dot(n_vector, e_vector)/(n*e))^2), dot(n_vector, e_vector)/(n*e));
f = atan2(sign(dot(r_vector, v_vector))*sqrt(1-(dot(r_vector,e_vector)/(r*e))^2),dot(r_vector,e_vector)/(r*e));
end