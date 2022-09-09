% Matty Hill
% ESS311
% MATLAB tutorial
% 01 / 04 / 22
% Section AA

% positional vector
positions = [0:40:400]

% velocity vector
vel = [0, 1.3, 3.1, 3.9, 4.1, 4.2, 4.3, 4.5, 4.8, 5.0, 5.1]

% Calculating movement total
tot_move = vel * 630

%plotting the vectors
landslide_plot = plot(positions, tot_move, "g-p")
title('Total Displacement of Landslide as a Function of Distance')
xlabel('Distance (meters)')
ylabel('Displacement (cm)')

% Saving the plot to a .jpeg file
saveas(landslide_plot, 'landslide_plot.jpeg')



