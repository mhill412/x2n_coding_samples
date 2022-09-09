function rms = calcRMS(observed,predicted);
%
% function rms = calcRMS(observed,predicted);
%
% Calculates the root mean square (RMS) misfit between an observed
% and predited set of values.  Low RMS suggests that the models fits
% the data well.
%
% Input
%	observed - a vector of observed values
%	predicted - a vector of values predicted by a model
%
% Output
%	rms - root mean square misfit
%
% Notes
%	Units should be the same for the observed and predicted.
%	Observed and predicted vectors must have the same length.
%

d=(observed - predicted);  % misfit
n=length(observed);
rms=  sqrt( sum(d.^2) / n );



 
