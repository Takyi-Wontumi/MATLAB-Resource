% This is a demonstration to plot the relationship between damping ratio and power
clc;close all; clear all;

% variables with  values
M_struct = 100;     % kg mass of machine
m_struct = 100;       % kg mass 
m_clothes = 1;      % kg mass of clothes
m_veh = 1;          % kg ~not needed
k_struct = 9.5E06;  % N/m
k_veh = 1;          % N/m

e_radius = 0.375;  % m radius
c_veh = 75;        % Ns/m ~not needed
c_struct = 0;      % Ns/m ~not needed

% Calculating damping ratios
struct_damping = c_struct / (2 * sqrt(k_struct * m_struct));
veh_damping = c_veh / (2 * sqrt(k_veh * m_veh));
total_damping = struct_damping + veh_damping;

% range of total damping ratios for plotting
total_dampingRatio = linspace(0, 1, 100);

% Natural frequency calculations
struct_freq = sqrt(k_struct / m_struct); % rad/s
veh_freq = sqrt(k_veh / m_veh);          % rad/s

% range of r_ratio (ω_veh / ω_struct) values
r_ratio = linspace(0, 2, 100);  % ratio range from 0 to 2

% Amplitude of vibration
Y_amplitude = 0.001; % meters, example value

% Calculate power normalized and resonant power
power_norm = (m_struct * total_damping * (Y_amplitude^2) * (r_ratio.^3) .* (struct_freq^3)) ./ ...
             ((1 - (r_ratio.^2)).^2 + (2 * total_damping * r_ratio).^2);

% Power when veh_freq equals struct_freq (r_ratio = 1, Resonance)
power_res = (m_struct * (Y_amplitude^2) * (struct_freq^3)) / (4 * total_damping);

% Plotting power vs ratio
figure;
plot(r_ratio, power_norm, 'b', 'LineWidth', 1.5);
hold on;
yline(power_res, 'r--', 'Power at Resonance', 'LineWidth', 1);
title('Power (W) vs Ratio (\omega_{veh} / \omega_{struct})');
xlabel('Ratio (\omega_{veh} / \omega_{struct})');
ylabel('Power (W)');
legend('Power vs. Ratio', 'Resonant Power');
grid on;

%{
% Plotting Damping Ratio vs ratio
figure;
plot(r_ratio, total_dampingRatio, 'g', 'LineWidth', 1.5); % Plotting damping ratio vs ratio
title('Damping Ratio (\zeta_{total}) vs Ratio (\omega_{veh} / \omega_{struct})');
xlabel('Ratio (\omega_{veh} / \omega_{struct})');
ylabel('Damping Ratio (\zeta_{total})');
grid on;
%}
%add eom that describe the system in x, v, a