% Load and run the Simulink model
params;
load('p53_mdm2_results.mat', 't', 'p53', 'Mdm2');
model_name = 'EKF_nonconstant_pi';

% Set simulation time
sim_time = 360;

% Run the Simulink model
simOut = sim(model_name, 'StopTime', num2str(sim_time));

% Extract results from Simulink logs
t_sim = simOut.get("tout");      % Time vector
p53_sim = simOut.get("p53_data");   % p53 concentration data    
Mdm2_sim = simOut.get("Mdm2_data");    % Mdm2 concentration data

% Interpolate Simulink data to match DDE-generated time vector
p53= interp1(t_sim, p53_sim, t, 'pchip'); % Shape-preserving cubic interp.
Mdm2 = interp1(t_sim, Mdm2_sim, t, 'pchip');

% Plot the results
figure(1);
plot(t, p53, 'b-', t, Mdm2, 'r-');
xlabel('Time');
ylabel('Concentration');
legend('p53', 'Mdm2');
title('Simulated oscillations of p53 and Mdm2, 360 seconds simulation time');

% Plot the phase plane trajectory
figure(2);
plot(p53, Mdm2, 'k-', 'LineWidth', 1.5);
hold on;

% Compute nullclines numerically
ks = 2; k1 = 2; tau = 3; K1 = 0.1;
k2 = 2; K2 = 1; dx = 1; dy = 1; n = 2;

p53_values = linspace(min(p53), max(p53), 100);
Mdm2_nullcline = k2 * p53_values.^n ./ (K2^n + p53_values.^n) / dy;
p53_nullcline = (ks - dx * p53_values) ./ (k1 * p53_values ./ (K1 + p53_values));

plot(p53_values, Mdm2_nullcline, 'b-', 'LineWidth', 1.5);
plot(p53_values, p53_nullcline, 'r-', 'LineWidth', 1.5);
xlabel('p53');
ylabel('Mdm2');
legend('Trajectory', 'Mdm2-nullcline', 'p53-nullcline');
title('Simulink Phase Plane Trajectory with Nullclines');
hold off;

% Save results
save('p53_mdm2_simulink_results.mat', 't', 'p53', 'Mdm2');