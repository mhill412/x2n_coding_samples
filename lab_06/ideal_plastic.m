function velocity=ideal_plastic(depth,average_velocity);
%
% velocity=ideal_plastic(depth,average_velocity); 
%
% Calculates velocity profile for an ideal plastic (plug flow). 
% All strain occurs at the base.
%
% Input
%	depth - depth below the surface [m], 0 is surface, positive down
%	average_velocity - average velocity  [m/yr]
%
% Output
%	velocity- velocity calculated at depth [m/yr]
%

velocity = average_velocity*ones(size(depth));






