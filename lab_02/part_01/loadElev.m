% loadElev.m
%
% Loads 10 meter DEM of Seattle into Matlab
%
% Output    x: row of x values [m]
%           y: row of y values [m]
%           z: 2D array of elevation values [m]
%

clear all;

nRows = 2742;
nCols = 1269;
dx = 10;  % spacing between points [meters]
x = 0:dx:(nCols-1)*dx;
y = 0:dx:(nRows-1)*dx;

% open data file
disp('Loading data...');
fid=fopen('seattleElev.bil','r');  
z=fread(fid,inf,'uint16')/10;      % dm -> m
z=reshape(z,nCols,nRows); % convert vector into an array
fclose(fid);

clear fid;
