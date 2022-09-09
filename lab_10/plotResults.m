function plotresults(a);

% plotresults(a)
%
% function for plotting results of groundwater model
%
%   a=0     plotting mode for runModel; plots 2nd figure with well data
%   a=1     plotting mode for practiceModel
%

close all;

% set this flag to zero when ready to fit well data
practice=a;  % =1 for practice mode, =0 for scenario mode

%-----------------------------

disp('Plotting results.');

p=load('head.dat');
wells=load('pumpwell.dat');

%must change spacing here, if grid resolution changed
x=0:0.125:3.;  % [km]
z=0:0.125:2.5;

fig1=figure(1);
fig1.Position=[100 200 700 400];
fig1.Name='Model Input and Outputs';

ax2=subplot(2,3,1);
cond=load('cond.dat');
%contourf(x,z,cond);
%pcolor(x,z,cond);
imagesc(x,z,cond);
colormap(ax2,jet);
colorbar('vert');
axis ij;
axis equal tight;
title({'Model input:','Conductivity (m/s)'});
ylabel('depth (km)');
xlabel('horizontal distance (km)');

ax3=subplot(2,3,4);
load bottom.dat;
load top.dat;
load right.dat;
load left.dat;
po=NaN*zeros(21,25);
po(1,:)=bottom;
po(end,:)=top;
po(:,1)=left;
po(:,end)=right;
cs=colormap('jet');
cs=[1 1 1;cs];
imagesc(x,z,flipud(po));
%pcolor(x,z,po);
colormap(ax3,cs);
colorbar('vert');
xlabel('horizontal distance (km)');
ylabel('depth (km)');
title({'Model Input:','Perimeter head levels (m)'});
axis ij;
axis equal tight;


%figure(2);
ax4=subplot('Position',[0.35 0.15 0.65 0.7]);
v = 0:2:100;
[c,h]=contour(x,z,p,v);
clabel(c,h);
xlabel('horizontal distance (km)');
ylabel('depth (km)');
%colormap(ax4,jet);
hold on;

[px,pz]=gradient(-p);
u = cond.*px;
v = cond.*pz;
quiver(x,z,u,v,2,'k-');

% tank position
a=[2.5, 2.65, 2.65, 2.5]; % horiz
b=[0.1, 0.1, 0.05, 0.05]; % depth
h=patch(a,b,3);

%Location of wells
plot(wells(1),wells(2),'r*');
text(wells(1),wells(2), ['Pumping' num2str(wells(3)) 'm^2/yr']);

[X,Y] = meshgrid(x,z);

% draw streamlines from leaking tank
s1 = streamline(X,Y,px,pz,2.5,0.05);
set(s1, 'LineWidth', 1.5);
s2 = streamline(X,Y,px,pz,2.65,0.1);
set(s2, 'LineWidth', 1.5);

text(2.3,-0.08,'Buried Tank');
text(0.3,-0.2,'River');

title({'Model Output:','Contours of Head (m)','Arrows of Darcy Velocity'});
grid on;
hold off;
axis ij;
axis equal;

%figure(4)
%surf(x,z,p) % example of 3-d map of the head
%title('Head')

if practice==1, disp('Practice mode.'); return; end

% compare model input & output to observations
fig2=figure(2);
fig2.Name='Comparison of Model and Observations';
fig2.Position=[120 220 800 400];

file={'well08_obs.dat'; 'well16_obs.dat'; 'well24_obs.dat'};  % edit list
for i=1:size(file);
    welldat=load(file{i});
    dp=size(welldat,1);  % row id corresponding to depth of well
    id=str2num(file{i}(5:6));  % column in cond corresponding to well
  
    subplot(1,3,i);
    plot(p(:,id),z,'b.-');
    title({file{i}(1:6),'Head Level (m)'});
    hold on;
    ii=~isnan(welldat(:,2));
    plot(welldat(ii,2),z(ii),'sr');
    hold off;
    axis ij;
    xlabel('Head Level [m]');
    legend('Model Prediction','Observation','Location','best');
    
    % add well locations to figure(1)
    figure(1)
    hold on
    text(x(id)-0.14,-0.2,file{i}(1:6)); % add well label
    plot([x(id) x(id)],[0 z(dp)],'k');  % draw well line
    hold off
    figure(2)
end
figure(1)
