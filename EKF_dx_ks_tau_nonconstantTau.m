%% EKF-based Parameter Estimation for the p53-Mdm2 Model (DDE version)
%
% Now estimating: d_x (decay of x(t)), k_s (production of x), and tau (delay).
% We keep k1, n, K1, k2, K2, d_y fixed.

clear; clc; close all;

%% Load synthetic or experimental data
%  Assume this loads t, p53, Mdm2 as column vectors of equal length
load("p53_mdm2_simulink_results.mat");  % Must contain: t, p53, Mdm2

N_time = length(t);
h      = t(2) - t(1);   % sampling timestep

%% --- Fixed Parameters ---
k1  = 2;       % feedback strength for y(t−tau)*x(t)/(K1 + x^n)
K1  = 0.1;     
n   = 2;       
k2  = 2;       % rate in second equation
K2  = 1;       
d_y = 1;       % decay rate of y(t)

%% --- EKF: Which Parameters We Estimate ---
% estimate [ d_x ; k_s ; tau ].
d_x_true = 1.0;
k_s_true = 2.0;
tau_true_sine = 3 + 0.1 * sin(2 * pi * (1/60) * t);

% Initial guesses
d_x_guess = 1.5;   
k_s_guess = 1.5;
tau_guess = 2.9; % in seconds 

%% Augmented State: [ x(t); y(t); d_x; k_s; tau ]
% Initialize the state from the first data point
x0 = p53(1);
y0 = Mdm2(1);

z_est = [ x0; 
          y0;
          d_x_guess;
          k_s_guess;
          tau_guess ];
nx = length(z_est);

%% Covariances
P = diag([0.01, 0.01, 0.01, 0.16, 0.2]);   % Initial guess of state covariance
Q = diag([0.0001, 0.0001, 0.002, 0.000705, 0.004]);  % Process noise 
R = diag([0.2,  0.2]);                % Measurement noise for [x; y]

%% For storing estimates over time
z_history = zeros(nx, N_time);
z_history(:,1) = z_est;

y_buffer = zeros(N_time,1);
y_buffer(1) = y0;

for k = 1:(N_time-1)
    
    %-----------------------------------------------------------
    % 1) RETRIEVE CURRENT STATE ESTIMATES
    %-----------------------------------------------------------
    xk   = z_est(1);
    yk   = z_est(2);
    dxk  = z_est(3);
    ksk  = z_est(4);
    tauk = z_est(5);
    
    %-----------------------------------------------------------
    % 2) GET DELAYED Y = y(t - tau)
    %-----------------------------------------------------------
    % Convert tau to an integer index shift (roughly):% Suppose are at step k -> have z_est(5) = tauk
    % and want y(t - tauk).
    
    delay_float = tauk / h;            % fractional offset
    delay_floor = floor(delay_float);  % integer part
    alpha       = delay_float - delay_floor;  % fraction in [0,1)
    
    idxA = k - delay_floor;       % first index
    idxB = k - delay_floor - 1;   % second index
    
    % clamping indices if they go out of range:
    if idxB < 1
        % Not enough history, clamp everything to y_buffer(1).
        y_delayed        = y_buffer(1);
        dy_delayed_dtau  = 0;
    else
        yA = y_buffer(idxA); 
        yB = y_buffer(idxB);
    
        % linear interpolation
        y_delayed = (1 - alpha)*yA + alpha*yB;
    
        % The partial derivative wrt tau => 1/h * (yB - yA)
        dy_delayed_dtau = (yB - yA)/h;
    end

    %-----------------------------------------------------------
    % 3) STATE PREDICTION (DISCRETE EULER)
    %    x'(t) = ks - k1*(y(t−tau)*x)/(K1 + x^n) - dx*x
    %    y'(t) = k2/(K2^n + x^n) - d_y*y
    %-----------------------------------------------------------
    % x'(t) = k_s - k1*( y_delayed*x )/(K1 + x^n ) - d_x*x
    frac_term = (k1 * y_delayed * xk)/(K1 + max(xk^n,1e-9));
    
    x_next = xk + h*( ksk - frac_term - dxk*xk );

    y_next = yk + h*( k2/(K2^n + max(xk^n, 1e-9)) ...
                      - d_y * yk );
                  
    % The parameters are treated as constant in the prediction:
    dx_next  = dxk; 
    ks_next  = ksk;
    tau_next = tauk;

    z_pred = [ x_next;
               y_next;
               dx_next;
               ks_next;
               tau_next ];
    
    %-----------------------------------------------------------
    % 4) JACOBIAN OF f w.r.t. z = [x, y, dx, ks, tau]
    %-----------------------------------------------------------
    denom_xn = (K1 + max(xk^n,1e-9));  % avoid zero
    frac     = (k1 * y_delayed * xk) / denom_xn;
    
    % partial derivatives wrt x:
    %   df_x/dx = 1 + h*d/dx( ... ) 
    %   where the inside is  - frac - dx*x + ks(does not depend on x) ...
    % but frac depends on x in multiple ways.
    
    % d frac/dx = k1 * y_delayed * [ (1 * denom_xn) - (xk * n*xk^(n-1)) ] / denom_xn^2
    %            = k1 * y_delayed * [ denom_xn - n*xk^n ] / denom_xn^2

    dFrac_dx = k1 * y_delayed * ...
               ( denom_xn - n*xk^n ) / (denom_xn^2 + 1e-18);
           
    F = eye(nx);  % start with identity

    % ---- row 1 is x_{k+1} partials ----
    % d x_next / d x
    F(1,1) = 1 + h*( - dFrac_dx - dxk );
    % d x_next / d y
    F(1,2) = 0;  % y appears only in y_delayed, which is treated as separate
    % d x_next / d dx
    F(1,3) = -h * xk;     % derivative of (- dx*x) wrt dx
    % d x_next / d ks
    F(1,4) = h;           % derivative of (ks) wrt ks
    % d x_next / d tau
    %   =  - h * [ (k1 * xk)/(K1 + xk^n ) ] * d[y_delayed]/d tau
    %   and  d[y_delayed]/d tau = - dy_delayed_dt
    % partial( x_{k+1} ) / partial( tau )
    %   = -h * [ (k1*xk)/(K1 + x^n) ] * partial( y_delayed )/ partial( tau )
    %   = -h * [ (k1*xk)/(K1 + x^n) ] * ( dy_delayed_dtau ).
    F(1,5) = -h * ( (k1*xk)/(K1 + xk^n) ) * dy_delayed_dtau;


    % ---- row 2 is y_{k+1} partials ----
    % y_next = yk + h*[ k2/(K2^n + x^n) - d_y*y ]
    % partial wrt x
    denom_xn2 = (K2^n + max(xk^n,1e-9));
    dTerm_dx  = - k2 * (n*xk^(n-1)) / denom_xn2^2;  % derivative of k2/(K2^n + x^n)
    F(2,1) = 1 + h*( dTerm_dx ); 
    % partial wrt y
    F(2,2) = 1 - h*d_y;
    % partial wrt dx
    F(2,3) = 0;  % dx doesn't appear in eqn for y
    % partial wrt ks
    F(2,4) = 0;  % ks doesn't appear in eqn for y
    % partial wrt tau
    %   y'(t) has no delayed y term (only x^n inside), so ∂/∂tau = 0
    F(2,5) = 0;

    % ---- row 3,4,5 for d_x, k_s, tau themselves are simply carried forward
    % d_x_next = d_x + 0 => partial derivative = 1 w.r.t d_x, and 0 w.r.t others
    % k_s_next = k_s + 0 => partial derivative = 1 w.r.t k_s
    % tau_next = tau + 0 => partial derivative = 1 w.r.t tau
    % The rest remain 0.
    %
    % So F(3,3) = 1, F(4,4) = 1, F(5,5) = 1, etc.  (Already set by eye(nx).)
    
    %-----------------------------------------------------------
    % 5) TIME UPDATE (Prediction) for covariance
    %-----------------------------------------------------------
    P = F * P * F' + Q;
    
    %-----------------------------------------------------------
    % 6) MEASUREMENT UPDATE
    %    We measure [ x; y ] directly => H = [I_2, 0_2x3]
    %-----------------------------------------------------------
    H = [1 0 0 0 0; 
         0 1 0 0 0];  
    z_meas = [ p53(k+1);
               Mdm2(k+1) ];
    
    y_tilde = z_meas - z_pred(1:2);   % measurement residual
    S       = H * P * H' + R;
    K_gain  = P * H' / S;
    
    z_est = z_pred + K_gain * y_tilde;
    P     = (eye(nx) - K_gain*H)*P;

    %-----------------------------------------------------------
    % 7) Store and also update the buffer
    %-----------------------------------------------------------
    z_history(:,k+1) = z_est;
    
    % store updated y for future delay usage
    y_buffer(k+1) = z_est(2);  
end

%% -- Plot the Results --
time = t;

% Extract the final estimates
x_est   = z_history(1,:);
y_est   = z_history(2,:);
dx_est  = z_history(3,:);
ks_est  = z_history(4,:);
tau_est = z_history(5,:);

figure; 
subplot(3,1,1)
plot(time,dx_est,'LineWidth',2); hold on
yline(d_x_true,'r--','LineWidth',1.5);
ylabel('d_x'); grid on; 
legend('Est','True','Location','Best');
title('Estimated d_x');

subplot(3,1,2)
plot(time,ks_est,'LineWidth',2); hold on
yline(k_s_true,'r--','LineWidth',1.5);
ylabel('k_s'); grid on; 
legend('Est','True','Location','Best');
title('Estimated k_s');

subplot(3,1,3)
plot(t, tau_est, 'LineWidth', 2); hold on
plot(t, tau_true_sine, 'r--', 'LineWidth', 1.5); 
xlabel('Time'); 
ylabel('\tau'); 
grid on; 
legend('Est', 'True', 'Location', 'Best');
title('Estimated \tau');

figure;
plot(time,x_est,'r--','LineWidth',2); hold on
plot(time,y_est,'b--','LineWidth',2);
plot(time,p53, 'r-','LineWidth',1.5)
plot(time,Mdm2,'b-','LineWidth',1.5)
xlabel('Time'); ylabel('Concentration');
legend('x_{est}','y_{est}','x_{meas}','y_{meas}','Location','Best');
title('States vs. Measurements');
grid on;

n_smooth = 250; 
n_smooth_ks = 350;

% Find the indices from N_time-n_smooth+1 to N_time:
N_time   = size(z_history, 2);
startIdx = max(1, N_time - n_smooth + 1);
startIdx_ks = max(1, N_time - n_smooth_ks +1);
idxRange = startIdx : N_time;
idxRange_ks = startIdx_ks : N_time;

% Compute the mean (smoothed) parameters:
d_x_smooth  = mean(z_history(3, idxRange));
k_s_smooth  = mean(z_history(4, idxRange_ks));
tau_smooth  = mean(z_history(5, idxRange));

disp(d_x_smooth);
disp(k_s_smooth);
disp(tau_smooth);

k1  = 2;     
n   = 2;      
K1  = 0.1;    
k2  = 2;      
K2  = 1;      
dy = 1;
ks = k_s_smooth;
dx = d_x_smooth;
tau = tau_smooth;

% Time span and initial conditions as before

model_name = 'PSO_p53_mdm2_simulation';
load_system(model_name);

% Set simulation time
sim_time = 360;

% Run the Simulink model
simOut = sim(model_name, 'StopTime', num2str(sim_time));

% Extract results from Simulink logs
t_est = simOut.get("tout");      % Time vector
p53_est = simOut.get("p53_data");   % p53 concentration data    
Mdm2_est = simOut.get("Mdm2_data");    % Mdm2 concentration data

% 5) Plot comparison
figure;
plot(t, p53, 'r-', 'LineWidth',1.5); hold on
plot(t, Mdm2, 'b-', 'LineWidth',1.5);
plot(t_est, p53_est, 'r--', 'LineWidth',2);
plot(t_est, Mdm2_est,'b--', 'LineWidth',2);
xlabel('Time'); ylabel('Concentration');
legend('p53 (true)','Mdm2 (true)','p53 (est)','Mdm2 (est)','Location','Best');
title('Comparison with EKF Parameter Estimates');
grid on;