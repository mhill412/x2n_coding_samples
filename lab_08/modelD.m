function T = modelD(z, b, crustal_heat_production,thermal_conductivity,surface_heat_flow,surface_temp)
% This function will model the geotherm with zero heat production through
% the crust.
% The inputs are (z, crustal thickness, crustal heat production, thermal conductivity, surface heat flow, surface temp)
T_coeff = (crustal_heat_production*b.^2)/thermal_conductivity;
T1 = T_coeff.*(1-exp(-z./b));
T2 = ((crustal_heat_production.*b./thermal_conductivity)+(surface_heat_flow./thermal_conductivity)).*z;
T = T1 - T2 + surface_temp;
plot( T, z);  % Plots the geotherm
end

