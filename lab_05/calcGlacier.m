% Matty Hill
% ESS311
% Glaciers
% 02 / 07 / 22
% Section AA

%% Part I
% Load the elevation data and plot the surface topography for the
% Athtabasca glacier
plotAtha;

% Using glacierProfile.m to make a topographic profile of the glacier.
% [dist_atha, z_atha] = glacierProfile(x_atha, y_atha, h_atha);

% Creating a plot of our profile
figure 
plot(dist_atha,z_atha,'.-'); 
xlabel('distance along profile (m)'); 
ylabel('elevation (m)'); 
%close all

% Saving the profile data to a .mat file
save points_atha dist_atha z_atha;

% Calculating the differential elevation of the profile
zdiff_atha = z_atha - z_atha(1);  

% Making a plot of the differential data
figure 
plot(dist_atha,zdiff_atha,'.-'); 
title('Differential Stress along Profile');
xlabel('distance along profile (m)'); 
ylabel('Differential Elevation (m)'); hold on;
%close all

% Plotting equation 3 from the lab sheet
ty = linspace(1e4, 1e5, 20);
rho = 900;
f = 0.7;
g = 9.81;
h_thickness = [];
figure; 
for i=1:length(ty)
    h_thickness = sqrt((2*ty(i)*dist_atha)/(f*rho*g));  
    plot(dist_atha, h_thickness);
    hold on
end
title('Differential Stress along Profile');
xlabel('distance along profile (m)'); 
ylabel('Thickness (m)'); 
plot(dist_atha,zdiff_atha,'.-'); hold off

% Calculating and plotting the misfit 
%close all;
figure;
for i=1:length(ty)
    misfitty = (zdiff_atha - sqrt((2*ty(i)*dist_atha)/(f*rho*g))).^2;
    sum_misfit = sum(misfitty);
    plot(ty(i), sum_misfit, '.');
    hold on
end
title('Misfit over Differing Yield Stresses')
xlabel('Yield Stress');
ylabel('Misfit')
hold off;

%% Part II
%close all;
% Loading in the data and plotting it forr Saskatchewan and Robson glaciers
%plotSask;
%plotRobs;

% Creating the glacial profiles
%[dist_sask,z_sask] = glacierProfile(x_sask,y_sask,h_sask);
%[dist_robs,z_robs] = glacierProfile(x_robs,y_robs,h_robs);

% Creating quick plots to confirm the profiles
%close all;
% Sask
figure 
plot(dist_sask,z_sask,'.-'); 
xlabel('distance along profile (m)'); 
ylabel('elevation (m)'); 
% Robs
figure 
plot(dist_robs,z_robs,'.-'); 
xlabel('distance along profile (m)'); 
ylabel('elevation (m)'); 

% Saving the profile data to .mat files
save points_sask dist_sask z_sask;
save points_robs dist_robs z_robs;

% Calculating surface slope dh/dx for both glaciers
slope_sask = calcSlope(dist_sask, z_sask);
slope_robs = calcSlope(dist_robs, z_robs);

% Calculating the height of the ice
h_ice_sask = (ty(8))./(f*rho*g.*slope_sask);
h_ice_robs = (ty(8))./(f*rho*g.*slope_robs);

% Calculating the bedrock height
z_bed_sask = z_sask - h_ice_sask;
z_bed_robs = z_robs - h_ice_robs;
%close all
% Creating the plot of observed ice surface elevation and bedrock topography
% Saskatchewan:
figure
plot(dist_sask, z_sask,'.-'); hold on
xlabel('distance along profile (m)'); 
ylabel('elevation (m)'); 
title('Saskatchewan Glacier Bedrock Elevation vs. Ice Elevation');
plot(dist_sask, z_bed_sask);
legend('Ice','Bedrock');
hold off;
% Robson:
figure
plot(dist_robs, z_robs,'.-'); hold on
xlabel('distance along profile (m)'); 
ylabel('elevation (m)'); 
title('Robson Glacier Bedrock Elevation vs. Ice Elevation');
legend();
plot(dist_robs, z_bed_robs);
legend('Ice','Bedrock');
hold off;
