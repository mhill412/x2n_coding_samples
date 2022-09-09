% loadData.m
%
% This script loads the file elev.dat into the Matlab workspace
%  and defines the longitude and latitude arrays
%
% Drew Stolar 1-14-5
%

clear all;

load elev.dat
elev=elev/1000;  % convert to km


lat = 20:5/60:50;
long = -135:5/60:-60;

% load state boundaries
stateLL=load('StateBoundaries.dat');

