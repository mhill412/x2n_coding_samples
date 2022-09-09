% Matty Hill
% ESS311
% Crustal Heat Flow
% 02 / 21 / 22
% Section AA

% Establishing Constants
[up_crust_craton, up_crust_basin] = deal(20000); % m
[low_crust_thick_craton, low_crust_thick_basin] = deal(20000, 10000); % m
[surf_heat_flow_craton, surf_heat_flow_basin] = deal(-45e-3, -90e-3); % W/m^2
[surf_temp_craton, surf_temp_basin] = deal(0); % degrees C
[therm_cond_craton, therm_cond_basin] = deal(2); % W/m degrees C
[heat_prod_granite_craton, heat_prod_granite_basin] = deal(1.8e-6); % W/m^3
solidus = [0 2000 4000 6000 8000 10000 15000 20000 30000 40000; 950 825 765 725 710 695 685 680 675 672];
%% Part I
% Craton
modelA(therm_cond_craton, surf_heat_flow_craton, surf_temp_craton); hold on;
modelB(heat_prod_granite_craton, therm_cond_craton, surf_heat_flow_craton, surf_temp_craton);
modelC(heat_prod_granite_craton, therm_cond_craton, surf_heat_flow_craton, surf_temp_craton);
modelD(heat_prod_granite_craton, therm_cond_craton, surf_heat_flow_craton, surf_temp_craton);
plot(solidus(2,:), solidus(1,:));
hold off;
legend('Model A', 'Model B', 'Model C', 'Model D', 'Solidus');
title('Crustal Heat Production Models for the Craton');
xlabel(['Temperature (' char(176) 'C)']);
ylabel('Depth (m)');
axis ij;

% Basin and Range
figure;
modelA(therm_cond_basin, surf_heat_flow_basin, surf_temp_basin); hold on;
modelB(heat_prod_granite_basin, therm_cond_basin, surf_heat_flow_basin, surf_temp_basin);
modelC(heat_prod_granite_basin, therm_cond_basin, surf_heat_flow_basin, surf_temp_basin);
modelD(heat_prod_granite_basin, therm_cond_basin, surf_heat_flow_basin, surf_temp_basin);
plot(solidus(2,:), solidus(1,:));
hold off;
legend('Model A', 'Model B', 'Model C', 'Model D', 'Solidus');
title('Crustal Heat Production Models for the Basin and Range');
xlabel(['Temperature (' char(176) 'C)']);
ylabel('Depth (m)');
axis ij;

%% Part II
solidus = [0 2000 4000 6000 8000 10000 15000 20000 30000 40000; 950 825 765 725 710 695 685 680 675 672];
% I have defined the solidus array at the top of the file, the line above
% is exemplary.
% In the above code, I have added the solidus to our plots.