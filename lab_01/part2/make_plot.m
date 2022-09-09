% Matty Hill
% ESS311
% MATLAB tutorial
% 01 / 04 / 22
% Section AA

% Loading the landslide data table
load landslides.dat

% Making the plot of landslide volume 
landslide_vol = semilogx(landslides(:,1), landslides(:,2), 'p')
grid on
title('Landslide Volume over Runout Distances')
xlabel('Landslide Volume (km^3)')
ylabel('Runout Distance (km)')

% Saving the plot as a jpeg
saveas(landslide_vol, 'landslide_vol.jpeg')

% Testing the calcL function
volume = landslides(:,1);
runout = calcL(volume);
