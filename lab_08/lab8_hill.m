% Matty Hill
% Lab 8 - Strength of the Lithosphere
% Section AA
% 02/28/2022

%% Part I
% Defining variables (a lot of them)
z = linspace(0, 100e+3, 1000); % Depth
[Ts_o, Ts_c, Ts_v] = deal(0,0,470); % Surface Temperature
[A_o, A_c, A_v] = deal(3e-6,1e-6,3e-7); % Heat Prodution coefficient
[qs_o, qs_c, qs_v] = deal(-70e-3,-58e-3,-15e-3); % Surface Heat Flow
[k_o, k_c, k_v] = deal(2,2.5,2); % Thermal Conductivity
[b_o, b_c, b_v] = deal(15000,40000,20000); % Crustal Thickness
[ESd_o, ESd_c, ESd_v] = deal(100); % Effective Stress Factor (dry)
[ESw_o, ESw_c, ESw_v] = deal(80,80,0); % Effective Stress Factor (wet)
[mu_o, mu_c, mu_v] = deal(0.6); % Coefficient of Friction
[rhoC_o, rhoC_c, rhoC_v] = deal(3300,2600,2700); % Crustal Density
[rhoM_o, rhoM_c, rhoM_v] = deal(3300); % Mantle Density
[Acrust_o, Acrust_c, Acrust_v] = deal(7e+4,1e-5,1e-9); % Power-Law fit parameter for Crust
[Acrust_wo, Acrust_wc, Acrust_wv] = deal(9e+5,5e-3,5e-7); % Power-Law fit parameter for Crust (wet conditions)
[Qcrust_o, Qcrust_c, Qcrust_v] = deal(5e+5,2e5,2e5); % Activation Energy for Creep (crust)
[nC_o, nC_c, nC_v] = deal(3); % Power Exponent for crust
[Amantle_o, Amantle_c, Amantle_v] = deal(7e+4); % Power-Law fit parameter for Mantle
[Amantle_wo, Amantle_wc, Amantle_wv] = deal(9e+5); % Power-Law fit parameter for Mantle
[Qmantle_o, Qmantle_c, Qmantle_v] = deal(5e+5); % Activation Energy for Creep (mantle)
[nM_o, nM_c, nM_v] = deal(3); % Power Exponent for mantle
[e_rate_o, e_rate_c, e_rate_v] = deal(1e-14); % Differential Strain Rate

% Using Model D to calculate the geotherm for the three lithospheres.
% Oceanic Lithosphere
T_o = modelD(z, b_o, A_o, k_o, qs_o, Ts_o); hold on;

% Continental Lithosphere
T_c = modelD(z, b_c, A_c, k_c, qs_c, Ts_c);

% Venusian Lithosphere
T_v = modelD(z, b_v, A_v, k_v, qs_v, Ts_v); hold off;

% Plotting Adjustments
axis ij;
title('Crustal Heat Production Models for the 3 Lithospheres');
xlabel(['Temperature (' char(176) 'C)']);
ylabel('Depth (m)');
legend('Oceanic Lithosphere', 'Continental Lithosphere', 'Venusian Lithosphere');

% Calculating the differential stress at failure from the surface to 100km
ds_o = strength(z, T_o, b_o, mu_o, ESd_o, rhoC_o, rhoM_o, Acrust_o, Qcrust_o, nC_o, Amantle_o, Qmantle_o, nM_o, e_rate_o);
ds_c = strength(z, T_c, b_c, mu_c, ESd_c, rhoC_c, rhoM_c, Acrust_c, Qcrust_c, nC_c, Amantle_c, Qmantle_c, nM_c, e_rate_c);
ds_v = strength(z, T_v, b_v, mu_v, ESd_v, rhoC_v, rhoM_v, Acrust_v, Qcrust_v, nC_v, Amantle_v, Qmantle_v, nM_v, e_rate_v);




%% Part II
% Calculating the differential stress at failure from the surface to 100 km
% for wet conditions
dsw_o = strength(z, T_o, b_o, mu_o, ESw_o, rhoC_o, rhoM_o, Acrust_wo, Qcrust_o, nC_o, Amantle_wo, Qmantle_o, nM_o, e_rate_o);
dsw_c = strength(z, T_c, b_c, mu_c, ESw_c, rhoC_c, rhoM_c, Acrust_wc, Qcrust_c, nC_c, Amantle_wc, Qmantle_c, nM_c, e_rate_c);
dsw_v = strength(z, T_v, b_v, mu_v, ESw_v, rhoC_v, rhoM_v, Acrust_wv, Qcrust_v, nC_v, Amantle_wv, Qmantle_v, nM_v, e_rate_v);


% Plot the lithospheric strength for dry and wet conditions as a
% function of depth

% Oceanic Lithosphere
figure;
plot(ds_o, z); hold on;
plot(dsw_o, z); hold off;
axis ij;
title('Oceanic Lithospheric Strength Models');
xlabel('Strength (MPa)');
ylabel('Depth (m)');
legend('Oceanic Lithosphere', 'Oceanic Lithosphere (wet)');
% Continental Lithosphere
figure;
plot(ds_c, z); hold on;
plot(dsw_c, z); hold off;
axis ij;
title('Continental Lithospheric Strength Models');
xlabel('Strength (MPa)');
ylabel('Depth (m)');
legend('Continental Lithosphere', 'Continental Lithosphere (wet)');
% Venusian Lithosphere
figure;
plot(ds_v, z);
% plot(dsw_v, z); hold off; % No such thing as a wet venusian lithosphere
axis ij;
title('Venusian Lithospheric Strength Models');
xlabel('Strength (MPa)');
ylabel('Depth (m)');
legend('Venusian Lithosphere');

% Calculating the total force
F_o = abs(sum(diff(1e6.*ds_o).*diff(z.*1e3)));
F_c = abs(sum(diff(1e6.*ds_c).*diff(z.*1e3)));
F_v = abs(sum(diff(1e6.*ds_v).*diff(z.*1e3)));
