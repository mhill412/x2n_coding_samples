% This script is where you parameterize your grounwater flow model.
% This script is for Part I of the lab.
% Modify the vectors below to define the head values along the model
% boundaries.

% There is no need to run relax.m or plotResults.m individually.
% They are automatically run at the end of this script.
clear all;  % clears workspace

%---------------------------------
% Step 1:
% Build up hydraulic head boundary conditions by defining 
% the head (meters) of the groundwater along the perimeter of the 
% model domain.  Each node in the model has a dimension of 0.125 km.
%
%            Top
%       L|---------|R
%       e|         |i
%       f|         |g
%       t|         |h
%        |---------|t
%           Bottom

% Replace the numbers below to build up left & right perimeter 
% (from bottom to top).
% For 0.125 km spacing, vector must be 21 elements long
right  = [2 1 0 -1 -2 -3 -4 -5 -6 -7 -8 -9 -10 -11 -12 -13 -14 -15 -16 -17 -18];
right = fliplr(right)+18;
left = [6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26];
left = left+18;

% Build up top & bottom perimeter (from left to right)
% For 0.125 km spacing, vector must be 25 elements long
top    = [26 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2]+18;
bottom = [6 5 4 3 2 1 0 -1 -2 -3 -4 -5 -6 -7 -8 -9 -10 -11 -12 -13 -14 -15 -16 -17 -18]+18;



%---------------------------------------------------------
% Nothing more to edit.  Ignore everything in the script below this line.
%----------------------------------------------------------

% Define the location and pumping rate for a well. Only for Part II.
% format: [x_location(km) z_location(km) pumping_rate(m^2/year)]
%in spatial coordinates, (0,0) is the upper left corner, positive downward
wells = [0 0 0];  % initialize well location

disp('Running practiceModel.m');
disp('Checking inputs.');
if size(top,2)~=25
      warning('Check length of the top vector. There should be 25 elements.'); return;
elseif size(bottom,2)~=25
    warning('Check length of the bottom vector. There should be 25 elements.'); return;
elseif size(left,2)~=21
    warning('Check length of the left vector. There should be 21 elements.'); return;
elseif size(right,2)~=21
    warning('Check length of the right vector. There should be 21 elements.'); return;
end

disp('Building conditions.');

% Build up cond by building each row (from left to right)
% hydraulic conductivity [m/s]
% 0.5 km spacing: 11 rows, 13 columns
% 0.25 km spacing: 21 rows, 25 columns: need to make more lines for this
cond = ones(21, 25);  % change this depending upon spacing
cond=5e-7*cond;  % uniform background value


% check that pumping well is in an allowed location
if wells(1)>0.25 && wells(1)<0.75
    warning('Pumping well is too close to the river.  Choose a different location.');
    return
elseif wells(1)>2.2 && wells(1)<2.8
    warning('Pumping well is within private property.  Choose a different location.');
elseif wells(2)<0
    warning('Pumping well z-coordinate is defined as positive downward.');
    return
end

% save as .dat files
save left.dat left -ascii
save right.dat right -ascii
save top.dat top -ascii
save bottom.dat bottom -ascii
save cond.dat cond -ascii
save pumpwell.dat wells -ascii

% clear all
relax  % Calculate finite difference solution of Laplace's equation for head

plotResults(1);  % plot results

disp('Done.');

