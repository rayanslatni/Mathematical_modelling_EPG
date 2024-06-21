%Find a function for findinding the fixed points of the system, in the
%simplified model (R,B)

% Define the default parameters
params_dict = struct(...
    'tau_B', 10, ...                  % timescale of BBB recovery [days]
    'tau_R', 10, ...                  % timescale of circuit remodeling [days]
    'k_IB', 0.1, ...                  % scaling parameter for the effect of neuroinflammation on BBB permeability [-]
    'k_BI', 1, ...                    % scaling parameter for the proinflammatory effect of BBB leakage [-]
    'k_IS', 2, ...                    % scaling parameter for the strength of seizure-promoting effects of neuroinflammation [-]
    'k_RS', 2, ...                    % scaling parameter for the strength of seizure-promoting effects of circuit remodeling [-]
    'k_ID', 8, ...                    % scaling parameter for the neurotoxic effect of overactivated glia [-]
    'k_BR', 1, ...                    % scaling parameter of BBB leakage on circuit remodeling [-]
    'k_DR', 0.0005, ...               % scaling parameter of neuronal loss on circuit remodeling [-]
    'K_SB', 0.875, ...                % scaling parameter for seizure burden on BBB integrity [-]
    'D_m', 1, ...                     % maximum possible extent of neuronal loss [-]
    'Theta', 0.25, ...                % Neurotoxicity threshold of overactivated glia [-]
    'IC', [-1,-1], ...           % initial conditions
    'D_const', 0,...               %Set D_cost as the critical value
    'IBDR_E_duration', [0, 2, 2, 0], ...   % Integration time step [days]
    'IBDR_E_amplitude', [0, 1.65, 1, 0], ... % Integration amplitude size
    'Complex_input', 'no', ...        % flag for complex simulation
    'amount_simulations', 4, ...      % amount of simulations
    'number_simulation', 1 ...       % number of simulations
);


% Add simulation-specific parameters
t_end = 5000; % simulation time
dt = 1;    % time step
IC = [-1 -1]; % initial conditions
t_vec = linspace(0, t_end, int16(t_end/dt + 1));
%a=-0.1;
%b=1;
%n=100000;
threshold_position = params_dict.Theta/params_dict.k_BI;
opt=odeset('RelTol',5e-6,'AbsTol',5e-4);
% ODE45 solution
[t, Y] = ode45(@(t, y) dIBDRdt_BBB_Rate_simplified(t, y, params_dict), t_vec, IC, opt);
B = Y(:,1);
R = Y(:,2);


%COMPUTE THE ISOCLINES
% up parameter values
B_values = linspace(-1, 2, 20);
R_values = linspace(-1, 2, 20);

% Create a grid of (B, R) values
[B_grid, R_grid] = meshgrid(B_values, R_values);

% Initialize arrays to store derivatives
dR_dt = zeros(size(R_grid));
dB_dt = zeros(size(B_grid));

% Evaluate derivatives at each point in the grid
for i = 1:numel(B_grid)
    x = [B_grid(i); R_grid(i)]; % Initial state vector
    dxdt = dIBDRdt_BBB_Rate_simplified(0, x, params_dict); % Compute derivatives
    dB_dt(i) = dxdt(1);
    dR_dt(i) = dxdt(2);
end

%COMPUTE THE EQUILIBRIUM
% % Define the system of equations
%equations = @(x) dIBDRdt_BBB_Rate_simplified(0, x, params_dict);
% 
% % Initial guess for the solution
%initial_guess = [0.3; 0.4]; % Adjust the initial guess as needed
% 
% % Solve the system of equations
%equilibrium_solution = fsolve(equations, initial_guess);
% 
%disp('Equilibrium Solution:');
%disp(['B: ', num2str(equilibrium_solution(1))]);
%disp(['R: ', num2str(equilibrium_solution(2))]);

figure(1);
hold off;
plot(B,R)
xlabel('BBB distruption');
ylabel('Remodelling');
title('B vs R in the simplified model');

%Plot the simplified model
figure(2);
hold on;
plot(t,R, 'Color', "#77AC30", "LineWidth",1) 
plot(t,B, 'Color', "#EDB120", "LineWidth",1)
grid on;
xlabel('Time (days)')
ylabel('Model variables')
title('EPG simplified model');
legend('B', 'R');

% Original color
originalColor =  [.5 .5 .5];
% Make it more transparent (e.g., set alpha to 0.5 for half transparency)
transparentColor = [originalColor(1:3), 0.5];

% Plot isoclines
figure(3);
hold on; 
quiver(B_grid, R_grid, dB_dt, dR_dt, 'AutoScaleFactor',1.9, 'Color',transparentColor);
contour(B_grid, R_grid, dB_dt,[2 0], "Color",[0.8,0.13, 0.6], "LineWidth",1.5); 
contour(B_grid, R_grid, dR_dt,[2 0], "Color",[0.2,0.13, 0.8], "LineWidth",1.5);
plot([0.41 0.41],ylim, 'r--', 'LineWidth', 1.5, 'DisplayName','Neurotoxicity threshold');  
%title('Isoclines and Direction Field in the (B, R): K_S_B = 0.875');
title('Isoclines and Direction Field in the (B, R): D_cost = 0.35 and K_D_R = 0.001 ');
xlabel('Extent of BBB distruption');
ylabel('Degree of circuit Remodelling');
legend('dB/dt Isoclines','dR/dt Isoclines' ,'Direction Field');
grid on;




%function for the simplified model
function dxdt = dIBDRdt_BBB_Rate_simplified(t, y, params_dict)
        B0 = y(1); R0 = y(2);
        %Define the f function
        f = params_dict.K_SB * (exp(params_dict.k_IS * (params_dict.k_BI * B0).^2 + params_dict.k_RS * R0) - 1) / ...
        (exp(params_dict.k_IS * (params_dict.k_BI * B0).^2 + params_dict.k_RS * R0) + 1);

        %definisci le ODE
        dBdt = 1/params_dict.tau_B*(-B0+params_dict.k_IB*params_dict.k_BI*B0+f);
        dRdt = 1./params_dict.tau_R*(-R0+params_dict.k_BR*B0+params_dict.k_DR*params_dict.D_const); %cambia il valore di D_cost
        %vector with solutions
        dxdt = [dBdt; dRdt]; % Ensure dxdt is a column vector

end





