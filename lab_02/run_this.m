% Matty Hill
% ESS311
% Slope Stability
% 01/11/22
% Section AA

%% load DEM data
loadElev;

%% Make maps of elevation and slopes
% Make 2 maps in one figure using the subplot function

% Plot elevation map
%subplot(1,2,1), plotElev;


% Calculate slope and plot slope map
calcSlope;
subplot(1,2,2), plotSlope;
caxis([0 30]);


%% Part 1
%% Calculate FS and make maps for both colluvium and glacial outwash
% Make 2 maps in one figure using the subplot function

% colluvium
FS_col = calcFS(s, 1, 1530, 5000, 0.45);
subplot(1,2,1), plotFS(x,y,z,FS_col);
title('Factor of Safety for Colluvium')
% glacial outwash
%FS_GO = calcFS(s, 1, 1796, 20000, 0.62);
%subplot(1,2,2), plotFS(x,y,z,FS_GO);
%title('Factor of Safety for Glacial Outwash')
%% Part 2
%% Calculate FS (pore pressure included) and make maps for colluvium material
% Make 2 maps in one figure using the subplot function

% colluvium, hw = 0.5
%FS_col_pp05 = calcFS_pp(s, 1, 1530, 5000, 0.45, 0.5);
%subplot(1,2,1), plotFS(x,y,z,FS_col_pp05);
%title('Factor of Safety With Pore Pressure = 0.5');
% colluvium, hw = 1
FS_col_pp1 = calcFS_pp(s, 1, 1530, 5000, 0.45, 1);
subplot(1,2,2), plotFS(x,y,z,FS_col_pp1);
title('Factor of Safety With Pore Pressure = 1');




