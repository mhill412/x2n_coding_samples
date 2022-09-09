% Marc Jaffrey and Zach Ploskey 2013
% This script computes the head map when given boundary conditions, a map 
%  of hydraulic conductivity, and a well location and pumping rate.  It produces 
%  an output file called head.dat, which contains the head values at the 
%  points where conductivity is defined.
%
% The only modification that you will ever have to make to this file is to
%  delx, which is the spacing of the head and conductivity grids.   It is 
%  set to a default value of 1000 m.  Be sure to change this if you use a
%  different spacing.
%

% Grid spacing in meters
delx = 125; 

%
% Do not change anything below this point!
%
toler = 0.00001;

%load the hydraulic conductivity data
load cond.dat
t = flipud(cond);
%t = cond;

% convert units from m/sec to m/yr so that it compares with pumping units
t=t*60*60*24*365;
% Add scale factor to get unts of finite diff calc correct- need to check.
% This use to be transmissivity, and thus multiplied by a thickness.
%t=t*1000;

%identify size
[nrow, ncol] = size(t);
ho = zeros(nrow, ncol);
q = zeros(nrow, ncol);

%load the boundary conditions
load top.dat
load bottom.dat
load left.dat
load right.dat

%Assign boundary condition values to old height matrix
[t1, t2] = size(top);
if t1 == 1
    ho(1, :) = bottom;
    ho(nrow, :) = top;
    ho(:, 1) = left;
    ho(:, ncol) = right;
else
    ho(1, :) = bottom';
    ho(nrow, :) = top';
    ho(:, 1) = left';
    ho(:, ncol) = right';
end

%enter one pumping well
%well location is given from upper left corner in x,z in kilometers

wells=load('pumpwell.dat');  % pumping units [m^2/yr]


ir = nrow - (wells(2) * 1000 + delx / 2) / delx + 1;  % depth (meters)
i = sscanf(num2str(ir), '%i');  % col index location for well
jr = (wells(1) * 1000 + delx / 2) / delx + 1;  % horiz distance (meters)
j = sscanf(num2str(jr), '%i');  % row index location for well
q(i, j) = q(i, j) - wells(3);   % insert point discharge rate at well


%main loop for finite difference solution
zr = zeros(nrow - 1, ncol - 1);
hl = zr; hr = zr; hu = zr; hd = zr;

diff = toler + 1;

td = t(1:end - 2, 2:end - 1);
tu = t(3:end, 2:end - 1);

tl = t(2:end - 1, 1:end - 2);
tr = t(2:end - 1, 3:end);

tc = t(2:end-1,2:end-1);

h = ho;
    
disp('Solving Laplaces equation...');

while (diff > toler)
    
    hd = h(1:end - 2, 2:end - 1);
    hu = h(3:end, 2:end - 1);    
    hl = h(2:end - 1, 1:end - 2);
    hr = h(2:end - 1, 3:end);
    
    % Flux solver form (K*Hx)x + (K*Hz)z=0       
    % -KDh=q (deals withthe well)
    
    h(2:end - 1, 2:end - 1) = (delx * q(2:end-1, 2:end-1) + (tc + tl) .* hl ...
        + (tc + tr) .* hr + (tc + tu) .* hu + (tc + td) .* hd) ...
        ./ (4 * tc + tl + tr + tu + td);

    diff = max(max(abs(h - ho)));

    ho = h;
end


h1 = flipud(h);

save head.dat -ascii h1
