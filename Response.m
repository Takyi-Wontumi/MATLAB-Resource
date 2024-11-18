
%m, M, r, c, k, F_init, forced_freq, x_cond, xdot_cond
clc;close all; clear all;
% Define parameters
m = 1; %kg
M = 100;%kg
r = 0.375; %m
c = 1000; %Ns/m
k = 9.5E06; %N/m
forced_freq = ((2*pi)/60)*(16E02); %ω converting 1600rpm to rad/s
disp(round(forced_freq,4));
F_init = m*r*forced_freq^2; %N
x_cond = 0; %m
xdot_cond = forced_freq*r; %m/s - converting angular freq to linear velocity v = ωr

% assigning response function
[displacement, velocity, acceleration, xResponse] = responseFunc(m, M, r, c, k, F_init, forced_freq, x_cond, xdot_cond);

% time range
time = linspace(0, 10, 1000);

% Evaluate the displacement, velocity, and acceleration over time
x_vals = displacement(time);
v_vals = velocity(time);
a_vals = acceleration(time);
totalResponse = xResponse;
disp(totalResponse);
fprintf('y(t) = %.3fsin(%.3ft)', displacement(0), forced_freq);

% Plotting results
figure;
subplot(3,1,1);
plot(time, x_vals);
title('Displacement Response');
xlabel('Time (s)');
ylabel('Displacement');

subplot(3,1,2);
plot(time, v_vals);
title('Velocity Response');
xlabel('Time (s)');
ylabel('Velocity');

subplot(3,1,3);
plot(time, a_vals);
title('Acceleration Response');
xlabel('Time (s)');
ylabel('Acceleration');
