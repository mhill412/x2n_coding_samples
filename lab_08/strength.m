function [diff_stress]=strength(z,T,dcrust,mu,eff_stress,den_crust,den_mantle,A_crust,Q_crust,n_crust,A_mantle,Q_mantle,n_mantle,e_rate);
%
% [diff_stress]=strength(z,T,dcrust,mu,eff_stress,den_crust,den_mantle,A_crust,Q_crust,n_crust,A_mantle,Q_mantle,n_mantle,e_rate);
%
%   Calculates the strength of the lithosphere for the brittle crust and viscous mantle.
%
%   Input Parameters [units]
%   ------------------------
%   z           depths, zero at surface, positive downward [meters]
%   T           temperature, defined at each depth z [Celsius]
%   dcrust	crustal thickness [m]
%   mu          coeficient of friction
%   eff_stress  effective stress factor [%].  100% indicates that 
%               the effective stress is equal to the normal stress (no pore 
%               pressure).  0% indicates an effective stress of zero.
%   den_crust	crustal density [kg/m^3]
%   den_mantle	mantle density  [kg/m^3]
%   A_crust     power-law fit parameter for the crust [MPa^(-3) s^(-1)]
%   Q_crust     activation energy for creep in the crust [J mol^(-1)]
%   n_crust     power exponnent in the crust
%   A_mantle	power-law fit parameter for the mantle [MPa^(-3) s^(-1)] 
%   Q_mantle    activation energy for creep in the mantle [J mol^(-1)]
%   n_mantle    power exponnent in the mantle 
%   e_rate      differential strain rate [strain/sec]
%
%
%   Output Parameters [units]
%   ------------------------
%   diff_stress     differential stress for failure	[MPa]
%
%
%   The frictional failure criteria assumes Byerlee's law:
%   tau_friction = mu * (sig_normal) + C 
%
%   The viscous failure criteria assumes a power-law viscous fluid:
%   e_rate = A * (sig_diff)^n * exp( -Q/R/T )
%

% print input parameters
disp(['crustal thickness [m]: ',num2str(dcrust)]);
disp(['coeficient of friction: ',num2str(mu)]);
disp(['effective stress factor [%]: ',num2str(eff_stress)]);
disp(['crustal density [kg/m^3]: ',num2str(den_crust)]);
disp(['mantle density [kg/m^3]: ',num2str(den_mantle)]);
disp(['power-law fit parameter for the crust [MPa^(-3) s^(-1)]: ',num2str(A_crust)]);
disp(['activation energy for creep in the crust [J mol^(-1)]: ',num2str(Q_crust)]);
disp(['power exponnent in the crust: ',num2str(n_crust)]);
disp(['power-law fit parameter for the mantle [MPa^(-3) s^(-1)]: ',num2str(A_mantle)]);
disp(['activation energy for creep in the mantle [J mol^(-1)]: ',num2str(Q_mantle)]);
disp(['power exponnent in the mantle: ',num2str(n_mantle)]);
disp(['differential strain rate [strain/sec]: ',num2str(e_rate)]);
disp(['------']);

% constants
R = 8.314 ; % ideal gas constant [J  K^(-1) mol^(-1) ]
g = 9.8; % gravity [m/s^2]

% calculate frictional strength
% eqt 5.36 from Stuwe (2002), assumes strike-slip faulting
sig_diff_fric_crust = 2* (6e7+mu*den_crust*g*z)*(eff_stress/100) / (sqrt(mu^2+1) + mu);  % [Pa]
sig_diff_fric_mantle = 2* (6e7+mu*den_mantle*g*z)*(eff_stress/100) / (sqrt(mu^2+1) + mu); % [Pa]
sig_diff_fric_crust = sig_diff_fric_crust/1e6;  % [MPa]
sig_diff_fric_mantle = sig_diff_fric_mantle/1e6; % [MPa]

% calculate power law strength
Tk=T+273.15; % Celcius to Kelvin
sig_diff_power_crust = exp(Q_crust/n_crust/R./Tk) .* (e_rate/A_crust).^(1/n_crust);  % [MPa]
sig_diff_power_mantle = exp(Q_mantle/n_mantle/R./Tk) .* (e_rate/A_mantle).^(1/n_mantle); %  [MPa]

icrust=find(z<=dcrust);  % find values in the crust
imantle=find(z>dcrust);  % find values in the mantle 
s=[sig_diff_fric_crust(icrust) sig_diff_fric_mantle(imantle) ; sig_diff_power_crust(icrust) sig_diff_power_mantle(imantle)]; 

% set output parameters
diff_stress=min(s,[],1);  % take the minimium strength at each depth [MPa]

end


