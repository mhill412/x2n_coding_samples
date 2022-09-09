% This script is where you parameterize your grounwater flow model.
% This is only for Part II of the lab.
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
left  = linspace(40,12,21);
%[7 9 11 13 15 17 17 18 20 22 24 26 28 30 30 30 30 30 30 30 30];
right = linspace(30,53,21);

% Build up top & bottom perimeter (from left to right)
% For 0.125 km spacing, vector must be 25 elements long
top    = [6 4 0 0 5 10 10 12 13 18 23 29 37 45 50 53 53 53 53 53 53 53 54 55 55];
bottom = linspace(30,30,25);


%---------------------------------
% Step 2: 
% Do not edit until Part II of lab. 
% Define the location and pumping rate for a well.
% format: [x_location(km) z_location(km) pumping_rate(m^2/year)]
%in spatial coordinates, (0,0) is the upper left corner, positive downward
%wells = [0 0 0];  
wells = [2.25 0.5 900];  % example parameters for pumping well

%------------------------------------------------------------------------
% Nothing more to edit.  Ignore everything in the script below this line.
%------------------------------------------------------------------------
disp('Running runModel.m');
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

cond=5e-7*cond;  % bckground
cond(1:9,:)=8e-5;  % upper layers
cond(3:6,3:12)=4e-4; % shallow  conductivity layer
cond(9:13,:)=6e-4;  % deep high conductivity layer


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

plotResults(0);  % plot results

disp('Done.');

