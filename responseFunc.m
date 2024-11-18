function [displacementResponse, velocityResponse, accelerationResponse, totalResponse] = responseFunc(m, M, r, c, k, F_init, forced_freq, x_cond, xdot_cond)
    % responseFuncSecond calculates the displacement, velocity, and acceleration responses
    % for a damped mass-spring system subject to forced vibration.
    % Inputs:
    %   m - mass of the system (kg)
    %   M - equivalent mass constant for response calculation
    %   r - response ratio parameter
    %   c - damping coefficient (Ns/m)
    %   k - spring constant (N/m)
    %   F_init - initial force (N)
    %   forced_freq - forcing frequency (rad/s)
    %   x_cond - initial displacement (m)
    %   xdot_cond - initial velocity (m/s)
    
    % Calculating extra variables needed for computations
    nat_freq = round(sqrt(k/M), 3);   % Natural frequency (rad/s)
    static_def = round(F_init/k, 3);  % Static deflection (m)
    damping_ratio =round( c / (2 * sqrt(k * M)), 3); % Damping ratio
    X_amplitude = round(((m * r * forced_freq^2) / sqrt((k - M * forced_freq^2)^2 + (2 * c * forced_freq)^2)), 3);
    

    % Calculate phase angle
    if c == 0
        phase = 0;
    else
        phase = atan((2 * damping_ratio * r) / (1 - r^2));
    end

    syms U V t  % U and V represent C1 and C2 respectively

    % Determine the total response based on the damping ratio
    if (c == 0) && (damping_ratio == 0)
        disp("Undamped system (ζ = 0)")
        X_amplitude = vpa((F_init / (k - m * forced_freq^2)), 3);
        x_particular = X_amplitude * sin(forced_freq * t);
        x_homogeneous = U * cos(nat_freq * t) + V * sin(nat_freq * t);

    elseif damping_ratio < 1
        disp("Underdamped system (ζ < 1)")
        damp_freq = nat_freq * sqrt(1 - damping_ratio^2); % Damped frequency
        x_homogeneous = exp(-damping_ratio * nat_freq * t) * ...
                        (U * cos(damp_freq * t) + V * sin(damp_freq * t));
        x_particular = X_amplitude * sin(forced_freq * t - phase);

    elseif damping_ratio > 1
        disp("Overdamped system (ζ > 1)")
        s1 = round(-damping_ratio * nat_freq - nat_freq * sqrt(damping_ratio^2 - 1), 3);
        s2 = round(-damping_ratio * nat_freq + nat_freq * sqrt(damping_ratio^2 - 1), 3);
        x_homogeneous = U * exp(s1 * t) + V * exp(s2 * t);
        x_particular = X_amplitude * sin(forced_freq * t - phase);

    elseif abs(damping_ratio - 1) < 1e-5  % Tolerance for critically damped
        disp("Critically damped system (ζ = 1)")
        x_homogeneous = (U + V * t) * exp(-nat_freq * t);
        x_particular = X_amplitude * sin(forced_freq * t - phase);

    else
        error('Unexpected damping ratio value.');
    end

    % Total displacement solution
    x_total = x_homogeneous + x_particular;

    % Initial conditions to solve for M and V
    eqn1 = subs(x_total, t, 0) == x_cond;
    eqn2 = subs(diff(x_total, t), t, 0) == xdot_cond;

    % Solve for constants M and V
    total_solution = solve([eqn1, eqn2], [U, V]);

    % Substitute M and V back into the total solution
    x_total = vpa(subs(x_total, [U, V], [total_solution.U, total_solution.V]),3);
    totalResponse = simplify(x_total);

    % Calculate velocity and acceleration by differentiating with respect to time
    velocity_total = diff(x_total, t);
    acceleration_total = diff(velocity_total, t);

    % Convert symbolic expressions into MATLAB functions with respect to t
    displacementResponse = matlabFunction(x_total, 'Vars', t);
    velocityResponse = matlabFunction(velocity_total, 'Vars', t);
    accelerationResponse = matlabFunction(acceleration_total, 'Vars', t);
    
    
    % Display variables with formatted output
    fprintf('--- Input Variables ---\n');
    fprintf('\n');
    fprintf('mass(m): %.3f kg |\t', m);
    fprintf('Mass(M): %.3f kg |\t', M);
    fprintf('Radius(r): %.3f m\t\n', r);
    fprintf('(ω_forced): %.3f rad/s |\t', forced_freq);
    fprintf('(k): %.3f N/m\n', k);
    fprintf('(c): %.3f Ns/m |\t', c);
    fprintf('Initial Force (F): %.3f N\n', F_init);
    fprintf('x_cond (x_o): %.3f m |\t', x_cond);
    fprintf('xdot_cond (xdot_o): %.3f m/s\n', xdot_cond);
    fprintf('\n');
    
    fprintf('--- Calculated Variables ---\n');
    fprintf('\n');
    fprintf('(ω_nat): %.3f rad/s |\t', nat_freq);
    fprintf('(δ): %.3f m | \t', static_def);
    fprintf('(ζ): %.3f\n', damping_ratio);
    fprintf('X: %.3f m\n', X_amplitude);
    fprintf('(ϕ): %.3f rad\n', phase);

    if damping_ratio < 1
        damp_freq = nat_freq * sqrt(1 - damping_ratio^2); % Recalculate if needed
        fprintf('(ω_d): %.3f rad/s\n', damp_freq);
    end

    fprintf('Initial Condition Constants:\n');
    fprintf('U (C1): %.3f |\t', vpa(total_solution.U, 3));
    fprintf('V (C2): %.3f\n', vpa(total_solution.V, 3));
    fprintf('---------------------------\n');

end
