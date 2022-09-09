% calcSlope
%
% Calculates the slope in the direction of steepest descent
%
% Input     z: 2D array of elevation values generated by loadElev [m]
%          dx: spacing of points in DEM (also generated by loadElev) [m]
%
% Output    s: slope [radians]
%

disp('Calculating slope of surface...');
[dzdx,dzdy] = gradient(z,dx);  % surface derivative
s = atan(sqrt(dzdx.^2+dzdy.^2));  % surface slope in radians
