% Define a model for EPG in the case of BBB disruption injury
% Using the deterministic rate for S(I,R)

% Define the default parameters
params_dict = struct(...
    'tau_I', 1, ...                   % timescale of neuroinflammatory reaction [days]
    'tau_B', 10, ...                  % timescale of BBB recovery [days]
    'tau_D', 10, ...                  % timescale of neuronal death process [days]
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
    'IC', [0, 0, 0, 0], ...           % initial conditions
    'IBDR_E_duration', [0, 7, 0, 0], ...   % Integration time step [days]
    'IBDR_E_amplitude', [0, 0.25, 0, 0], ... % Integration amplitude size
    'Complex_input', 'no', ...        % flag for complex simulation
    'amount_simulations', 4, ...      % amount of simulations
    'number_simulation', 1 ...       % number of simulations
);

% Add simulation-specific parameters
t_end = 500; % simulation time
dt = 0.1;    % time step
IC = [0 0 0 0]; % initial conditions
t_vec = linspace(0, t_end, int16(t_end/dt + 1));
opt=odeset('RelTol',5e-6,'AbsTol',5e-4);
% ODE45 solution
[t, Y] = ode45(@(t, y) dIBDRdt_BBB_Rate(t, y, params_dict), t_vec, IC, opt);

% Extract the variables
I_vec = Y(:, 1); % Inflammation
B_vec = Y(:, 2); % BBB disruption
D_vec = Y(:, 3); % Neuron loss
R_vec = Y(:, 4); % Remodeling

%Dictionary with results
results_dict = struct('t_vec', t_vec, 'I_vec', I_vec, 'B_vec', B_vec,'D_vec', D_vec, 'R_vec', R_vec);
%Save the results
filename = strcat("C:\Users\rayan\OneDrive\Documenti\MATLAB\Mathematical Modeling\EPG");
save(filename, 'params_dict', 'results_dict');

%PLOT
figure(1);
hold on;
plot(t,I_vec, 'Color', "#D95319", "LineWidth",1) 
plot(t,B_vec, 'Color', "#EDB120", "LineWidth",1)
plot(t,D_vec, 'Color', "#7E2F8E", "LineWidth",1)
plot(t,R_vec, 'Color', "#77AC30", "LineWidth",1)
grid on;
xlabel('Time')
ylabel('Model variables')
title('EPG model');
legend('I', 'B', 'D', 'R')

%PLOT R vs B
figure(2);
hold off;
plot(B_vec,R_vec, 'Color', "#D95319", "LineWidth",1) 
xlabel('Extent of BBB distruption')
ylabel('Degree of circuit Remodelling')
title('B vs R');



% Define the function for the BBB model
function dxdt = dIBDRdt_BBB_Rate(t, y, params_dict)
    % Extract the initial condition
    I0 = y(1); B0 = y(2); D0 = y(3); R0 = y(4);
    inp = zeros(size(params_dict.IC));

    for ww = 1:length(inp)
        % calculating external input at this timestep.
        % Input may be simple (from time step 0 to T_off)
        if t <= params_dict.IBDR_E_duration(ww)
            inp(ww) = params_dict.IBDR_E_amplitude(ww);
        end
    end

    S = SeizureBurden(I0, R0, params_dict);
    dIdt = (1 / params_dict.tau_I) * (-I0 + params_dict.k_BI * B0 + inp(1));
    dBdt = (1 / params_dict.tau_B) * (-B0 + params_dict.k_IB * I0 + S + inp(2));
    dDdt = (1 / params_dict.tau_D) * ((1 - D0 / params_dict.D_m) * params_dict.k_ID * ...
        max(0, I0 - params_dict.Theta) + inp(3));
    dRdt = (1 / params_dict.tau_R) * (-R0 + params_dict.k_BR * B0 + params_dict.k_DR * D0 + inp(4));
    dxdt = [dIdt; dBdt; dDdt; dRdt]; % Ensure dxdt is a column vector
end

% Define the seizure burden function (Eq.4) depending on neuroinflammation and circuit remodeling
function SB = SeizureBurden(I, R, params_dict)
    SB = params_dict.K_SB * (exp(params_dict.k_IS * I * I + params_dict.k_RS * R) - 1) / ...
        (exp(params_dict.k_IS * I * I + params_dict.k_RS * R) + 1);
end

% Define the fun_dbdt function
function result_dbdt = fun_dbdt(x, params_dict)
    result_dbdt = -x + params_dict.k_IB * params_dict.k_BI * x + params_dict.K_SB * ...
        (exp(params_dict.k_IS * params_dict.k_BI * params_dict.k_BI * x * x + ...
        params_dict.k_RS * params_dict.k_BR * x + params_dict.k_RS * params_dict.k_DR * ...
        params_dict.D_const) - 1) / (exp(params_dict.k_IS * params_dict.k_BI * ...
        params_dict.k_BI * x * x + params_dict.k_RS * params_dict.k_BR * x + ...
        params_dict.k_RS * params_dict.k_DR * params_dict.D_const) + 1);
end


% a=-0.1;
% b=1
% n=100000
% 
% %Define a function for finding the equilibrium
% function output = fixedpointsfinder(result_dbdt, a, b, n)
%     output = cell(0, 2);
%     x = linspace(a, b, n+1);
% 
%     for ii = 2:(length(x)-1)
%         if fun(x(ii)) <= 0
%             if result_dbdt(x(ii-1)) > 0 && result_dbdt(x(ii+1)) < 0
%                 output{end+1, 1} = x(ii);
%                 output{end, 2} = 'Stable';
%             elseif result_dbdt(x(ii-1)) < 0 && result_dbdt(x(ii+1)) > 0
%                 output{end+1, 1} = x(ii);
%                 output{end, 2} = 'Unstable';
%             end
%         elseif result_dbdt(x(ii)) >= 0
%             if result_dbdt(x(ii-1)) < 0 && result_dbdt(x(ii+1)) < 0
%                 output{end+1, 1} = x(ii);
%                 output{end, 2} = 'Semistable';
%             elseif result_dbdt(x(ii-1)) > 0 && result_dbdt(x(ii+1)) > 0
%                 output{end+1, 1} = x(ii);
%                 output{end, 2} = 'Semistable';
%             end
%         end
%     end
% end

