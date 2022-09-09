% loadSeismic.m
%
%
% This script loads the file seismicThickness.dat into the MATLAB workspace
%
% The data is an estimate of the crustal thickness from the travel
% time of seismic waves.
%


seisT=load('seismicThickness.dat');
seisT=seisT/1000;  % convert to km
% seisT= flipud(seisT);  % uncomment for section AA

% load state boundaries
stateLL=load('StateBoundaries.dat');

