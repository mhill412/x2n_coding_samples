% Matty Hill
% ESS311
% Constitutive Laws
% 02 / 14 / 22
% Section AA

% Loading your data
data = load('vel_profile.dat');

% Generating a plot using errorbar
errorbar_plot = errorbar(data(:,2), data(:,1), data(:,3));
set(gca, 'YDir','reverse');
title('Salmon Glacier Uncertainty');
xlabel('Velocity (m/yr)');
ylabel('Depth (m)');

% Establishing constants
rho_ice = 900; % kg/m^3
surf_slope = 2; % degrees
thickness = 500; % m
velocity = data(:,2);
depth = data(:,1);
ice_visc = 1e13; % Pa s
power_law = 3;
yield_stress = 5e4; % Pa

% Calculating each of the velocity profiles
vel_newtonian = newtonian(depth, ice_visc, rho_ice, surf_slope, thickness);
vel_power_no_slip = power_no_slip(depth, power_law, thickness, max(velocity));
vel_power_slip = power_slip(depth, power_law, rho_ice, surf_slope, thickness, min(velocity));
vel_ideal_plastic = ideal_plastic(depth, mean(velocity));
vel_bingham_plastic = bingham_plastic(depth, ice_visc, yield_stress, rho_ice, surf_slope, thickness);

% Generating a plot of all velocity profiles for the inital parameters
plot(vel_newtonian, depth); hold on
plot(vel_power_no_slip, depth); hold on
plot(vel_power_slip, depth); hold on
plot(vel_ideal_plastic, depth); hold on
plot(vel_bingham_plastic, depth); hold on
plot(velocity, depth); hold off
set(gca, 'YDir','reverse');
title('Velocity Profiles for Each Prediction and Observations');
xlabel('Velocity (m/yr)');
ylabel('Depth (m)');
legend('Newtonian', 'Power No Slip', 'Power Slip', 'Ideal Plastic', 'Bingham Plastic', 'Observations', 'Location', 'northwest');

% Adjusting the parameters for optimized fit
ice_visc_newt = 1.5e13; % Pa s
ice_visc_bing = 7e12; % Pa s 7e12
power_law_no_slip = 3.8;
power_law_slip = 3.032;
yield_stress = 5e4; % Pa

% Calculating any new velocity profiles for the optimized fit
vel_newtonian = newtonian(depth, ice_visc_newt, rho_ice, surf_slope, thickness);
vel_power_no_slip = power_no_slip(depth, power_law_no_slip, thickness, max(velocity));
vel_power_slip = power_slip(depth, power_law_slip, rho_ice, surf_slope, thickness, min(velocity));
vel_bingham_plastic = bingham_plastic(depth, ice_visc_bing, yield_stress, rho_ice, surf_slope, thickness);

% Generating a plot of all velocity profiles for the adjusted parameters
figure;
plot(vel_newtonian, depth); hold on
plot(vel_power_no_slip, depth); hold on
plot(vel_power_slip, depth); hold on
plot(vel_ideal_plastic, depth); hold on
plot(vel_bingham_plastic, depth); hold on
plot(velocity, depth); hold off
set(gca, 'YDir','reverse');
title('Velocity Profiles for Each Prediction and Observations');
xlabel('Velocity (m/yr)');
ylabel('Depth (m)');
legend('Newtonian', 'Power No Slip', 'Power Slip', 'Ideal Plastic', 'Bingham Plastic', 'Observations', 'Location', 'northwest');

% Calculating the RMS misfit for each predicted curve
newt_rms = calcRMS(velocity, vel_newtonian);
power_no_rms = calcRMS(velocity, vel_power_no_slip);
power_rms = calcRMS(velocity, vel_power_slip);
ideal_rms = calcRMS(velocity, vel_ideal_plastic);
bing_rms = calcRMS(velocity, vel_bingham_plastic);



