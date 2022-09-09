%%% Part II
%% 5.
% Loading in the data set
d = load('corona_list.txt');
diameter = d(:,1); % first column is the diameter [km]
corona_diameter = 475; % km

% Adding the estimated diameter to the dataset
diameter(end+1)= corona_diameter; 

% Converting from diameter (km) to radius (m)
radius = diameter.*500; 

% Calculating the mean and standard deviation of the observations
r_avg = mean(radius);
r_std = std(radius);

% Plotting a histogram of the radii
figure('Name','hist_radii','NumberTitle','off');
histogram(radius, 20);
xlabel('Radius (m)');
ylabel('Number of Instances');
title('Histogram of Radii of Venusian Coronas');

%% 6.
% Calculating the vertical velocity 
t_diff = 200; % degrees Celcius
rho = 4500; % kg/m^3
eta = 1e+21; % Pa s
alpha = 3e-5; % 1/C
g = 8.87; % m/s^2
v = ((alpha*rho*g.*(radius.^2)*t_diff)./(3*eta)).*3.1536e+7; % m/year

% Creating a histogram of velocities
figure('Name','hist_velocity','NumberTitle','off');
histogram(v, 100);
xlabel('Vertical Velocity (m/year)');
ylabel('Number of Instances');
title('Histogram of Vertical Velocities for Venusian Coronas');

% Calculating the mean and standard deviation of the velocities
v_avg = mean(v);
v_std = std(v);
