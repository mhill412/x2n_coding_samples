function T = modelB(crustal_heat_production,thermal_conductivity,surface_heat_flow,surface_temp)
% This function will model the geotherm with zero heat production through
% the crust.
% The inputs are (crustal heat production, thermal conductivity, surface heat flow, surface temp)
z = linspace(0, 40000, 41000);
T = -(crustal_heat_production./(2.*thermal_conductivity)).*(z.^2) -(surface_heat_flow./thermal_conductivity).*z + surface_temp;
plot( T, z);  % Plots the geotherm
end

