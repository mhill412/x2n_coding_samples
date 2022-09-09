function []=plotMohrCircle(sigma_x,sigma_y,tau_xy)
%
% function []=plotMohrCircle(sigma_x,sigma_y,tau_xy)
%
% Plot the 2D Mohr Cicle of Stress.
%
% sigma_x - max normal stress acting on a plane [MPa]
% sigma_y - min normal stress acting on a conjugate plane [MPa]
% tauxy - shear stress acting on both planes [MPa]
%
% It is assumed that sigma_x > sigma_y.

if sigma_x < sigma_y
	disp('Error: sigma_x must be greater than sigma_y.');
	return
end

gridsize=1000;  % controls how many points are plotted along circle
phi=linspace(0,pi,gridsize);
sigma_mohr=(sigma_x+sigma_y)/2+(sigma_x-sigma_y)/2*cos(2*phi)+...
    tau_xy*sin(2*phi);
tau_mohr=-(sigma_x-sigma_y)/2*sin(2*phi)+...
    tau_xy*cos(2*phi);
sigma_1=(sigma_x+sigma_y)/2-sqrt(((sigma_x-sigma_y)/2)^2+tau_xy^2);
sigma_2=(sigma_x+sigma_y)/2+sqrt(((sigma_x-sigma_y)/2)^2+tau_xy^2);
tau_1=sqrt(((sigma_x-sigma_y)/2)^2+tau_xy^2);
tau_2=-tau_1;
center_circle=(sigma_x+sigma_y)/2;

plot(sigma_mohr,tau_mohr);
ax=axis;
fac=0.1; % scale factor for expanding window of figure
dx=fac*(ax(2)-ax(1));
dy=fac*(ax(4)-ax(3));
axis([ax(1)-dx ax(2)+dx ax(3)-dy ax(4)+dy]);
grid on;
axis equal;
xlabel('Normal Stress, MPa');
ylabel('Shear Stress, MPa');

hold on;
plot(sigma_1,0,'r*',sigma_2,0,'r*',...
    center_circle,tau_1,'ro',...
    sigma_x,tau_xy,'k.',sigma_y,-tau_xy,'k.');
% hold off




