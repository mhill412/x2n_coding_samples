function [slope]=calcSlope(distance,elevation)
%
% function [slope]=calcSlope(distance,elevation)
%
% Calculates surface slope
%
% Input:
%	distance - distance along linear profile [m]
%	elevation - surface elevation along linear profile [m]
%
% Ouput:
%	slope - surface slope along linear profile
%
%
% Notes:
% - slope calculated at midpoint, smoothed, and interpolated at 
%   original points along profile.
% - it is assumed that points are in the proper sequence along profile
% - it is also assumed that elevation increases along the profile
%

% mid-point between postings
distanceMid = distance(1:end-1) + diff(distance)/2;

% initial estimate of slope using differential operator
s=diff(elevation)./diff(distance);
s(1)=1;  % seed the first slope to be a high number

% apply moving average
s_smooth=movingmean(s',15);

% interpolate slope at original points along the profile
slope=interp1(distanceMid,s_smooth,distance);



