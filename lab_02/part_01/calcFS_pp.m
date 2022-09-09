function FoS_pp = calcFS_pp(slope,depth,density,cohesion,f,hw)
     
% function FoS = calcFS(s,depth,density,cohesion,f)
%
% Calculate the factor of safety (FoS) for a dry slope
%
% Input 
%       slope: 2D array of slope values calculated by calcSlope [radians]
%       depth: Depth to the failure plane [m]
%     density: Density of the material [kg/m3]
%    cohesion: Cohesion of the material [Pa]
%           f: Coef of Friction of the material [unitless]
%
% Output   FoS: Factor of safety
%


%%%%%%   ADD EQUATION HERE  %%%%%
% Tip: To avoid errors from typographical errors, consider calculating the 
% numerator and denomerator on separate lines.
    FoSnum = cohesion + (f.*(density.*9.8.*depth.*cos(slope).*cos(slope)-(997.*9.8.*hw)));
    FoSden = density.*9.8.*depth.*cos(slope).*sin(slope);
    FoS_pp = FoSnum./FoSden;
end