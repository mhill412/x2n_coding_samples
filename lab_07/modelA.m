function T = modelA(thermal_conductivity,surface_heat_flow,surface_temp)
% This function will model the geotherm with zero heat production through
% the crust.
% The inputs are (thermal conductivity, surface heat flow, surface temp)
z = linspace(0, 40000, 41000);
T = -(surface_heat_flow./thermal_conductivity).*z + surface_temp;
plot( T, z);  % Plots the geotherm
end

