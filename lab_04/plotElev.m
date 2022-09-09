% plotElev.m
%
% This script plots the elevation data
%
%

% plot elevation data
imagesc(long,lat,elev); hold on;

% add contour of coastline
%contour(long,lat,elev,[0 0],'k'); hold off;

% add state boundaries
plot(stateLL(:,2),stateLL(:,1),'k'); 
hold off;

% label the axes
xlabel('Longitude [^o E]');
ylabel('Latitude [^o N]');

% reset the axes
axis xy
axis([-135 -60 20 50]);

% add a colobar
colorbar;

% set aspect ratio of map
daspect([1 .67 1])

% add title
title('Elevation above sea level, [km]');


