function velocity=newtonian(depth,viscosity,density,slope,thickness);
%
% velocity=newtonian(depth,viscosity,density,slope,thickness);
%
% Calculates velocity profile for a Newtonian fluid
%
% Input
%	depth- depth below the surface [m], 0 is surface, positive down
%	viscosity- dynamic viscosity of material [Pa s]
%	density- density of material [kg/m^3]	
%	slope- surface slope [degrees], 0 is horiztonal
%	thickness - total thickness of body [m]
%
% Output
%	velocity- velocity calculated at depth [m/yr]
%

% constants
g=9.8;  % m/s^2, gravity
spyr=60*60*24*365; % s/yr,  seconds per year

velocity = (density*g*sind(slope)/2/viscosity)*( 2*depth*thickness-depth.^2 )*spyr;

velocity=flipud(velocity);




