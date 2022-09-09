function T = modelC(crustal_heat_production,thermal_conductivity,surface_heat_flow,surface_temp)
% This function will model the geotherm with zero heat production through
% the crust.
% The inputs are (crustal heat production, thermal conductivity, surface heat flow, surface temp)
z_uppercrust = linspace(0,20000,21000); % upper crust depth values [meters] 
z_lowercrust = linspace(20000,40000,21000); % lower crust depth values
z = [ z_uppercrust  z_lowercrust ];  % concatenate depth vectors 
tuc = -(crustal_heat_production./(2.*thermal_conductivity)).*(21000.^2) -(surface_heat_flow./thermal_conductivity).*21000 + surface_temp;
dtuc = -(crustal_heat_production/thermal_conductivity).*21000 - (surface_heat_flow./thermal_conductivity);
T_uppercrust = -(crustal_heat_production./(2.*thermal_conductivity)).*(z_uppercrust.^2) -(surface_heat_flow./thermal_conductivity).*z_uppercrust + surface_temp;
T_lowercrust = dtuc.*(z_lowercrust-21000)+tuc;
T = [ T_uppercrust  T_lowercrust ];  % concatenate temperature vectors 
plot( T, z);  % Plots the geotherm
end


