function velocity=power_no_slip(depth,power,thickness,top_velocity);
%
% velocity=power_no_slip(depth,power,thickness,top_velocity);
%
% Calculates velocity profile for a Power Law fluid with no slip allowed 
% at bottom boundary.
%
% Input
%	depth - depth below the surface [m], 0 is surface, positive down
%	power - n exponent in power law where (strain rate) ~ (stress)^n 
%	thickness - total thickness of body [m]
%	top_velocity - velocity at the surface [m/yr]
%
% Output
%	velocity - velocity calculated at depth [m/yr]
%

% constants
xp=power+1;

velocity=top_velocity - top_velocity*(1-depth/thickness).^xp;

velocity=flipud(velocity);







