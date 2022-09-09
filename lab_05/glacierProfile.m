function [dist,z] = glacierProfile(x_grid,y_grid,h_grid)

% function [dist,z] = glacierProfile(x_grid,y_grid,h_grid)
%
% Makes topo profiles from glacier elevation datasets.
% Converts x-y coordinates to linear distance along profile.
%
% Input:
% x_grid - Vector of x-axis values for figure (km)
% y_grid - Vector of y-axis values for figure (km)
% h_grid - 2D matrix of height values for figure (m)
%
% Output:
% dist - Output linear distance along transect (m)
% z - Output heights along profile (m)
%
% Notes:
% - Takes user input from active figure. 
% - Requires (x,y,z) of figure matches input (x_grid,y_grid,h_grid). 
% - Press enter to stop choosing points and exit function.
%


% Select points along profile
disp('Start picking points on active figure. Press enter when done.')
[xi,yi]=ginput();


% Get variables ready...
dl = 0.04; %Point spacing for interpolation along profile line (km).
z = []; %Initialize height variable
dist = []; %Initialize distance variable


%Loop to interpolate values segment by segment
for i = 1:(length(xi)-1)
    % Pick segment start and end points
    a = [xi(i) xi(i+1)];
    b = [yi(i) yi(i+1)];
    
    % Set up offset distance for segment start point
    if  i==1
        offset=0;
    else
        offset=max(dist);
    end
    
    % Build x and y coords along segment using dx spacing
    L = sqrt((a(1)-a(2))^2 + (b(1)-b(2))^2);    
    d_seg = offset+(0:dl:L);
    dx = abs(a(1)-a(2))/(L/dl);
    dy = abs(b(1)-b(2))/(L/dl);
    x_seg = min(a)+(0:dx:(max(a)-min(a)));
    y_seg = min(b)+(0:dy:(max(b)-min(b)));
    
    % Flip order of vectors if needed
    if a(1) > a(2)
        x_seg = fliplr(x_seg);
    end
    
    if b(1) > b(2)
        y_seg = fliplr(y_seg);
    end
    
    % Interpolate z data along segment
    z_seg = griddata(x_grid,y_grid,h_grid,x_seg,y_seg);
    
    % Append to full profile dataset
    dist = [dist d_seg];
    z = [z z_seg];
end

% remove duplicate points along profile
bad=[];
for i=2:size(dist,2)
    if dist(i)==dist(i-1)
        bad=[bad; i];
    end
end

dist(bad)=[];
z(bad)=[];

%Convert to meters
dist=dist*1000;
