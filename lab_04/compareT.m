% Matty Hill
% ESS311
% Isostasy
% 01 / 31 / 22
% Section AA

%% Part I
% Loading in the data for the lab and plotting the map of the United States
loadElev;
plotElev;

% Calculating the thickness of the crustal root, H
H = (2500/800).*(elev+9.5);

% Calculating the total crust thickness, T
T = (elev + 9.5) + 5 + H;

% New figure with a map of total crustal thickness
figure;
imagesc(long,lat,T); hold on;

% add state boundaries
plot(stateLL(:,2),stateLL(:,1),'k'); 

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
title('Crustal Thickness, [km]');

% Assign coordinates to a variable 
% Interior Lowlands
x_IL = [-97.8937   
   -93.7598  
   -89.2323   
   -84.1142   
   -85.6890 
   -82.3425 
   -92.3819
   -93.1693 
   -88.6417 
   -90.6102];
y_IL = [47.6614
   44.2323
   43.3091
   42.5177
   39.7480
   40.4075
   39.8799
   42.3858
   40.8031
   45.9469];
% Basin and Range 
x_BR = [-120.3346
   -118.1693
   -115.6102
   -114.8228
   -118.5630
   -116.2008
   -116.3976
   -115.0197
   -113.4449
   -118.1693];
y_BR = [42.1220
   41.4626
   41.0669
   38.8248
   39.7480
   39.4843
   37.3740
   36.5827
   35.6594
   36.7146];
% Colorado Plateau
x_CP = [-110.0984
   -107.9331
   -105.3740
   -106.5551
   -107.3425
   -108.3268
   -110.6890
   -109.1142
   -109.9016
   -108.3268];
y_CP = [38.4291
   38.8248
   37.9016
   36.3189
   37.5059
   36.5827
   36.0551
   35.5276
   34.4724
   34.2087];
% Gulf Coastal Plain 
x_GCP = [-97.5000
   -93.9567
   -90.6102
   -89.0354
   -85.2953
   -81.9488
   -81.3583
   -79.3898
   -77.4213
   -81.9488];
y_GCP = [29.3287
   30.3839
   29.9882
   31.0433
   30.9114
   28.6693
   26.6909
   33.9449
   35.5276
   31.8346];
% Southern Rocky Mountains
x_SRM = [-105.3740
   -105.5709
   -106.7520
   -106.1614
   -104.5866
   -105.9646
   -105.3740
   -106.1614
   -106.7520
   -104.1929];
y_SRM = [41.9902
   40.6713
   40.4075
   40.0118
   38.9567
   38.5610
   37.6378
   37.3740
   38.2972
   38.6929];

x = [x_IL, x_BR, x_CP, x_GCP, x_SRM];
y = [y_IL, y_BR, y_CP, y_GCP, y_SRM];

% Plotting our points on the crustal thickness map
plot(x,y,'.');
hold off
% Using the interp2 function to interpolate the crustal thickness at each
% point
T_interp = interp2(long, lat, T, x, y);

% Now we will use mean.m and std.m to find each at every province. 
avg = mean(T_interp, 1);
sd = std(T_interp, 1);
names = ["Interior Lowlands", "Basin and Range", "Colorado Plateau", "Gulf Coastal Plain", "Southern Rocky Mts."];
varnames = ["Province", "Mean", "Standard Deviation"];
% Create the table
topo_table = table(names(:), avg(:), sd(:), 'VariableNames', varnames)

% Closing all previous tables
close all;

%% Part II
% Load in and plot the seismic data
loadSeismic;
plotSeismic;

% Interpolating the seismic data and calculating the means and standard
% deviations
T_seis_interp = interp2(long, lat, seisT, x, y);
avg_seis = mean(T_seis_interp, 1);
sd_seis = std(T_seis_interp, 1);

% Creating the table for the seismic data calculations
seis_table = table(names(:), avg_seis(:), sd_seis(:), 'VariableNames', varnames)

% Making a new plot that compares the data from my tables
figure;
errorbarxy(avg, avg_seis, sd, sd_seis, {'o','r','r'} ); hold on;
plot([30 65], [30 65], 'b--'); hold off;
xlim([42 60]);

% Adding text to label which province is which
text(avg+0.25, avg_seis-0.5, names);

% Adding axes names and a title
title('Comparison of Crustal Thickness with Error as Generated from Elevation and Seismic Data');
xlabel('Crustal Thickness from Topographical Data');
ylabel('Crustal Thickness from Seismic Data');

% Closing all figures again
close all;