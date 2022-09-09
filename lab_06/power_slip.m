function velocity=power_slip(depth,power,density,slope,thickness,bottom_velocity);
%
% velocity=power_slip(depth,power,density,slope,thickness,bottom_velocity);
%
% Calculates velocity profile for a Power Law fluid with slip allowed 
% at bottom boundary.
%
% Input
%	depth - depth below the surface [m], 0 is surface, positive down
%	power - n exponent in power law where (strain rate) ~ (stress)^n 
%	density - density of material [kg/m^3]	
%	slope - surface slope [degrees], 0 is horiztonal
%	thickness - total thickness of body [m]
%	bottom_velocity - bottom boundary condition [m/yr]
%
% Output
%	velocity- velocity calculated at depth [m/yr]
%

% constants
g=9.8;  % m/s^2, gravity
spyr=60*60*24*365; % s/yr,  seconds per year
a=2e-24; % 1/(Pa^3 s)
xp=power+1;


velocity=bottom_velocity + ...
 (a*spyr/xp) * ((density*g*sind(slope))^power) * ( thickness^xp - (thickness-depth).^xp );

velocity=flipud(velocity);






