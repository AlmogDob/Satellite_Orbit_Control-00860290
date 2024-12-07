function [f] = calculate_drag_force(r_vec, v_vec)
global_variabels;

r = norm(r_vec);

h = r - R_earth;

if (h > 1000) 
    f = 0;
else
    table = readtable("density_table.xlsx");
    
    for index = 1:length(table.base_height) 
        if (table.base_height(index) < h && table.ceil_height(index) > h)
            break;
        end
    end
    rho_b = table.base_density(index);
    h_b   = table.base_height(index);
    H     = table.scale_height(index);
    
    rho = rho_b*exp(-(h-h_b)/(H));
    
    K_D = A_drag*C_D/m;
    f = -0.5*rho*norm(v_vec)^2*K_D*v_vec/norm(v_vec);
end
end