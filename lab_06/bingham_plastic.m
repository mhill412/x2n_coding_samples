function velocity=bingham_plastic(depth,viscosity,yield_stress,density,slope,thickness);
%
% velocity=bingham_plastic(depth,viscosity,yield_stress,density,slope,thickness);
%
% Calculates velocity profile for a Bingham plastic (plastic-viscous body) 
%
% Input
%	depth- depth below the surface [m], 0 is surface, positive down
%	viscosity - dynamic viscosity of material [Pa s]
%	yield_stress - plastic yield stress of material [Pa]
%	density - density of material [kg/m^3]	
%	slope - surface slope [degrees], 0 is horiztonal
%	thickness - total thickness of body [m]
%
% Output
%	velocity- velocity calculated at depth [m/yr]
%

% constants
g=9.8;  % m/s^2, gravity
spyr=60*60*24*365; % s/yr,  seconds per year

% calculate velocity for lower viscous portion
velocity = (spyr*density*g*sind(slope)/viscosity)*(depth*thickness - 0.5*depth.^2 )  - spyr*yield_stress*depth/viscosity;

% calculate for upper plastic portion,
% find depths below yield stress where no strain will ocurr
i=find( (thickness-depth) <= yield_stress/density/g/sind(slope) );  % plug flow
velocity(i)=spyr/viscosity * (0.5*density*g*sind(slope)*thickness^2 + (yield_stress^2)/2/density/g/sind(slope) - yield_stress*thickness);

velocity=flipud(velocity);





