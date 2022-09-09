%
% plotSask
%
%
% Plots topographic map Sasketchewan glacier
%

display('Loading data...');
load data_saskatchewan

% Plotting map of Sasketchewan Glacier
display('Plotting Sasketchewan Glacier...');
figure
imagesc(x_sask,y_sask,h_sask);
hold on
[c,p]=contour(x_sask,y_sask,h_sask,[0:100:4000],'k');
clabel(c,p);
hold off
caxis([1500 3500]);
colormap(cmap);
colorbar('horiz');
axis equal tight
xlabel('distance east [km]');
ylabel('distance south [km]');
title('Sasketchewan Glacier, elevation [m]');
text(7,3.8,'Sasketchewan Glacier','Rotation',40,'color','white');
