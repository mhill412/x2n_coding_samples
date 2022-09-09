% plotSlope
%
% Plots slope of DEM
%
% Required variables:
%           x: row of x values generated by loadElev [m]
%           y: row of y values generated by loadElev [m]
%           z: 2D array of elevation values generated by loadElev [m]
%           s: 2D array of slope values generated by calcSlope [radians]
%

% Convert from radians to degrees
sDegrees = s*180/pi;

% Make image
disp('Plotting surface slope...');
%figure(2)
imagesc(x/1000,y/1000,sDegrees'); hold on;

% Set the axis properties
axis equal; axis([min(x) max(x) min(y) max(y)]/1000);

% Add sea level contour
contour(x/1000,y/1000,z',[1 1],'w');

xlabel('[km]'); ylabel('[km]');

% Make a colorbar;
caxis([0 65]);  % set the min,max of the color spectrum
cb = colorbar; 

% Set the title of the colorbar
cbt = get(cb,'title');
set(cbt,'string','Slope [^o]');