% plotSeismic.m
%
% This script plots the crustal thicknes imaged by seismic tomography
%
%

figure
% plot elevation data 
imagesc(long,lat,seisT); hold on; 

% add contour of coastline 
%contour(long,lat,elev,[0 0],'k'); hold off; 

% add state boundaries
plot(stateLL(:,2),stateLL(:,1),'k'); hold off;

% label the axes 
xlabel('Longitude [^o E]'); 
ylabel('Latitude [^o N]'); 
% reset the axes 
axis xy;  % comment out line for section AA
axis([-135 -60 20 50]); 
 
% add a colobar 
colorbar; 
caxis([0 65]);
% set aspect ratio of map 
daspect([1 .67 1]) 
% add title 
title('Crustal Thickness from Seismic Tomography, [km]'); 





