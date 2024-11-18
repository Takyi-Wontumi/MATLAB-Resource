%{
% Define the time range 
syms t
t_s = linspace(0, 10, 1000);  % Adjust the time range as needed

% Define the first function (example: linear function f1(t) = 2*t + 1)
f1 = 0.3547*exp(-2.249*t)*cos(15.65*t + 1.428) + 0.071*cos(7*pi + 0.4);  % Replace with your first function

% Define the second function (example: linear function f2(t) = 2*t + 1)
f2 = (71*cos(7*pi*t + 2/5))/1000 + exp(-(2533230379510295*t)/1125899906842624)*(sin((4405219871792607*t)/281474976710656)*((2498090418307072*pi*sin(2/5))/78664640567725125 - (35971871389046189*cos(2/5))/3524175897434085600 + 226699919596230977/704835179486817120) - cos((4405219871792607*t)/281474976710656)*((71*cos(2/5))/1000 - 3/200));  % Replace with your second function

% Evaluate the functions over the time range
f1_vals = double(subs(f1, t, t_s));  % Evaluate f1 for each value in t_s
f2_vals = double(subs(f2, t, t_s));  % Evaluate f2 for each value in t_s

% Plot both functions on the same graph
figure;
plot(t_s, f1_vals, 'b', 'LineWidth', 2);  % f1 in blue
hold on;
plot(t_s, f2_vals, 'r--', 'LineWidth', 2);  % f2 in red dashed
hold off;

% Add labels and legend
xlabel('Time (t)');
ylabel('Function Value');
title('Comparison of Two Functions');
legend('Function 1', 'Function 2');

% Check for equality by subtracting one from the other
difference = f1_vals - f2_vals;

% Display a message based on the difference
if all(abs(difference) < 1e-01)  % A small tolerance to account for numerical errors
    disp('The functions are equal.');
else
    disp('The functions are NOT equal.');
end


%m, M, r, c, k, F_init, forced_freq, x_cond, xdot_cond
clc;close all; clear all;
% Define parameters
m = 100; %kg
M = 100;%kg
r = 0.375; %m
c = 50; %Ns/m
k = 9.5E06; %N/m
forced_freq = ((2*pi)/60)*(16E02); %ω converting 1600rpm to rad/s
disp(round(forced_freq,4));
F_init = m*r*forced_freq^2; %N
x_cond = 0; %m
xdot_cond = forced_freq*r; %m/s - converting angular freq to linear velocity v = ωr

% assigning response function
[displacement, velocity, acceleration, xResponse] = responseFuncSecond(m, M, r, c, k, F_init, forced_freq, x_cond, xdot_cond);

% time range
time = linspace(0, 10, 1000);

% Evaluate the displacement, velocity, and acceleration over time
x_vals = displacement(time);
v_vals = velocity(time);
a_vals = acceleration(time);
totalResponse = xResponse;
disp(totalResponse);

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

%}

%{
% Define the time range 
syms t
t_s = linspace(0, 10, 1000);  % Adjust the time range as needed

% Define the first function (example: linear function f1(t) = 2*t + 1)
f1 = 0.3547*exp(-2.249*t)*cos(15.65*t + 1.428) + 0.071*cos(7*pi + 0.4);  % Replace with your first function

% Define the second function (example: linear function f2(t) = 2*t + 1)
f2 = (71*cos(7*pi*t + 2/5))/1000 + exp(-(2533230379510295*t)/1125899906842624)*(sin((4405219871792607*t)/281474976710656)*((2498090418307072*pi*sin(2/5))/78664640567725125 - (35971871389046189*cos(2/5))/3524175897434085600 + 226699919596230977/704835179486817120) - cos((4405219871792607*t)/281474976710656)*((71*cos(2/5))/1000 - 3/200));  % Replace with your second function

% Evaluate the functions over the time range
f1_vals = double(subs(f1, t, t_s));  % Evaluate f1 for each value in t_s
f2_vals = double(subs(f2, t, t_s));  % Evaluate f2 for each value in t_s

% Plot both functions on the same graph
figure;
plot(t_s, f1_vals, 'b', 'LineWidth', 2);  % f1 in blue
hold on;
plot(t_s, f2_vals, 'r--', 'LineWidth', 2);  % f2 in red dashed
hold off;

% Add labels and legend
xlabel('Time (t)');
ylabel('Function Value');
title('Comparison of Two Functions');
legend('Function 1', 'Function 2');

% Check for equality by subtracting one from the other
difference = f1_vals - f2_vals;

% Display a message based on the difference
if all(abs(difference) < 1e-01)  % A small tolerance to account for numerical errors
    disp('The functions are equal.');
else
    disp('The functions are NOT equal.');
end
%}

%m, M, r, c, k, F_init, forced_freq, x_cond, xdot_cond
clc;close all; clear all;
% Define parameters
m = 150; %kg
M = 150;%kg
r = 0; %m
c = 2000; %Ns/m
k = 25E03; %N/m
forced_freq = 20; %ω converting 1600rpm to rad/s
disp(round(forced_freq,4));
F_init = 100; %N
x_cond = 0.1; %m
xdot_cond = 1.5; %m/s - converting angular freq to linear velocity v = ωr

% assigning response function
[displacement, velocity, acceleration, xResponse] = responseFuncSecond(m, M, r, c, k, F_init, forced_freq, x_cond, xdot_cond);

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
