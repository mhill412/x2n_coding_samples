%
% plotAtha
%
%
% Plots topographic map for Athabasca Glacier
%


display('Loading data...');
load data_athabasca

% Plotting map of Athabasca Glacier
display('Plotting Athabasca Glacier...');
figure
imagesc(x_atha,y_atha,h_atha);
hold on
[c,p]=contour(x_atha,y_atha,h_atha,[0:100:4000],'k');
clabel(c,p);
hold off
caxis([1500 3500]);
colormap(cmap);
colorbar;
axis equal tight
xlabel('distance east [km]');
ylabel('distance south [km]');
title('Athabasca Glacier, elevation [m]');
text(2,3.5,'Athabasca Glacier','Rotation',45,'color','white');





