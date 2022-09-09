%
% plotRobs
%
%
% Plots topographic map for Robson glacier
%

display('Loading data...');
load data_robson

% Plotting map of Robson Glacier
display('Plotting Robson Glacier...');
figure
imagesc(x_robs,y_robs,h_robs);
hold on
[c,p]=contour(x_robs,y_robs,h_robs,[0:100:4000],'k');
clabel(c,p);
hold off
caxis([1500 3500]);
colormap(cmap);
colorbar
axis equal tight
xlabel('distance east [km]');
ylabel('distance south [km]');
title('Robson Glacier, elevation [m]');
text(3,1.5,'Robson Glacier','Rotation',-60,'color','white');




