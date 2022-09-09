% Matty Hill
% ESS311
% Fractures and State of Stress
% 01 / 18 / 22
% Section AA

%% First we are going to be dealing with the granite dataset
% Loading data in
data_granite=load('friction_experiment_granite.dat');

% Plotting the normal stress and shear stress as a scatterplot 
plot([data_granite(:,1)], [data_granite(:,2)], 'k.');
% Turning the grid on
grid on;
% Axes labels and title
xlabel('Normal Stress (MPa)'), ylabel('Shear Stress (MPa)');
title('Normal Stress and Shear Stress of Granite and Fault Gouge');

% Utilizing polyfit function to fit a line to lab data
p_granite = polyfit( [data_granite(:,1)], [data_granite(:,2)], 1); 

% Utilizing the polyval function to predict the points
y_predicted_granite = polyval( p_granite, [data_granite(:,1)]);

% Using the plot command to add a line to our figure
hold on
plot([data_granite(:,1)], y_predicted_granite, 'r');
text(53, 30, '\uparrow Granite')

%% Now the above steps will be repeated, but with the gouge dataset instead
% Loading data in
data_gouge=load('friction_experiment_fault_gouge.dat');

% Plotting the normal stress and shear stress as a scatterplot 
hold on
plot([data_gouge(:,1)], [data_gouge(:,2)], 'k.');

% Utilizing polyfit function to fit a line to lab data
p_gouge = polyfit( [data_gouge(:,1)], [data_gouge(:,2)], 1); 

% Utilizing the polyval function to predict the points
y_predicted_gouge = polyval( p_gouge, [data_gouge(:,1)]);

% Using the plot command to add a line to our figure
hold on
plot([data_gouge(:,1)], y_predicted_gouge, 'b');
text(77, 13.5, '\uparrow Fault Gouge')
hold off;

%% Now we are on part II
% First, we will plot our mohr circles
figure; hold on;
plotMohrCircle(86, 20, 0); hold on;
plotMohrCircle(129, 45, 0); hold on;
plotMohrCircle(180, 75, 0); hold on;
plotMohrCircle(222, 100, 0); hold on;
plot([0 180], [20 68.3]);
hold off;

% Now, we will make another figure with both the Mohr-Coulomb envelope and
% Beyerlee's Law
figure; hold on
mc_plot = plot([0 180], [4.0355 95.8175]); M1 = "Mohr-Coulomb";
beyer_plot = plot([0 180], [20 68.3]); M2 = "Beyerlee's"; hold on;
plotMohrCircle(72, 60, 0); hold on;
plotMohrCircle(110, 60, 0); hold on;
plotMohrCircle(154, 60, 0); hold on;
plotMohrCircle(184, 60, 0);
legend([mc_plot,beyer_plot], [M1, M2]);
title("Mohr-Coulomb Failure Envelope and Beyerlee's Law for Fractured Granite")
hold off;