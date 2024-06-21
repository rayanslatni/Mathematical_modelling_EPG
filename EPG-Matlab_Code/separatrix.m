%Role of Neuronal loss (D) in the dynamics
%plot 7b

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
D_vis_vect = [0, 0.3, 0.4103035, 1]; %4 different value for D



FP = cell(length(D_vis_vect), 1); % Initialize cell array to store fixed points

for ii = 1:length(D_vis_vect)
    D_const = D_vis_vect(ii);

    B_max = 1.1;
    R_max = 1.1;
    B_min = -0.1;
    R_min = -0.1;

    [R, B] = meshgrid(linspace(R_min, R_max, 5000), linspace(B_min, B_max, 5000));

    exp_term = exp(params_dict.k_IS * (params_dict.k_BI * B).^2 + params_dict.k_RS * R);
    f = params_dict.K_SB * (exp_term - 1) ./ (exp_term + 1);
    U = 1 / params_dict.tau_B * (-B + params_dict.k_IB * params_dict.k_BI * B + f); % equation of B
    V = 1 / params_dict.tau_R * (-R + params_dict.k_BR * B + params_dict.k_DR * D_const); % linear R
    FP{ii} = fixedpointsfinder(@fun_dbdt, -0.1, 1, 100000, params_dict, D_const);
    % velocity = sqrt(U .* U + V .* V);
end 


%FP = fixedpointsfinder(@fun_dbdt, -0.1, 1, 100000, params_dict); % Check for steady states on given interval with given discretization step
marker_sizer = 2.5; % Replace 2 with the appropriate value



%----
for cc = 1:length(FP)
    % Access the fixed points for the current D_const
    current_FP = FP{cc};

    % Check the condition using the first element of the fixed point
    if current_FP{1} >= params_dict.Theta
        col = 'grey'; % EPG state
    else
        col = 'black'; % healthy state
    

    mrksize = 10;
    if strcmp(current_FP{2}, 'Unstable') % white one
        plot(ax0, current_FP{1}, current_FP{1} + params_dict.k_DR * D_vis_vect(cc), 'o', 'Color', col, 'MarkerSize', mrksize/marker_sizer);
        plot(ax0, current_FP{1}, current_FP{1} + params_dict.k_DR * D_vis_vect(cc), 'o', 'Color', 'white', 'MarkerSize', (mrksize-3)/marker_sizer);
    elseif strcmp(current_FP{2}, 'Semistable')
        plot(ax0, current_FP{1}, current_FP{1} + params_dict.k_DR * D_vis_vect(cc), 'o', 'Color', col, 'MarkerSize', mrksize/marker_sizer);
        plot(ax0, current_FP{1}, current_FP{1} + params_dict.k_DR * D_vis_vect(cc), 'o', 'Color', 'white', 'MarkerSize', (mrksize-3)/marker_sizer);
        plot(ax0, current_FP{1}, current_FP{1} + params_dict.k_DR * D_vis_vect(cc), 'x', 'Color', col, 'MarkerSize', (mrksize-3)/marker_sizer);
    elseif strcmp(current_FP{2}, 'Stable')
        plot(ax0, current_FP{1}, current_FP{1} + params_dict.k_DR * D_vis_vect(cc), 'o', 'Color', col, 'MarkerSize', mrksize/marker_sizer);
    end
    end
end


% PLOTS
B_max = 0.1;
R_max = 0.1;
B_min = -0.01;
R_min = -0.01;

FP = cell(length(D_vis_vect), 1); % Initialize cell array to store fixed points

% Loop over D_vis_vect
for ii = 1:length(D_vis_vect)
    D_const = D_vis_vect(ii);

    % Create meshgrid
    [R, B] = meshgrid(linspace(R_min, R_max, 5000), linspace(B_min, B_max, 5000));

    % Calculate differential equations
    exp_term = exp(params_dict.k_IS * (params_dict.k_BI * B).^2 + params_dict.k_RS * R);
    f = params_dict.K_SB * (exp_term - 1) ./ (exp_term + 1);
    U = 1 / params_dict.tau_B * (-B + params_dict.k_IB * params_dict.k_BI * B + f); % equation of B
    V = 1 / params_dict.tau_R * (-R + params_dict.k_BR * B + params_dict.k_DR * D_const); % linear R
    FP{ii} = fixedpointsfinder(@fun_dbdt, -0.1, 1, 100000, params_dict, D_const);

    % Find fixed points
    %FP = fixedpointsfinder(fun_dbdt, -0.1, 1, 100000); % You need to define fun_dbdt

    % Plotting
    figure;
    strm = quiver(B, R, U, V, 'linewidth', 1, 'color', [0.75, 0.75, 0.75], 'density', [0.5, 1.5]);
    hold on;

    for cc = 1:length(FP)
        if FP(cc, 1) >= params_dict.Theta
            col = 'grey';
        else
            col = 'black';
        end

        mrksize = 10;
        if strcmp(FP(cc, 2), 'Unstable') || strcmp(FP(cc, 2), 'Semistable')
            plot(FP(cc, 1), FP(cc, 1) + params_dict.k_DR * D_const, 'o', 'color', col, 'markersize', mrksize);
            plot(FP(cc, 1), FP(cc, 1) + params_dict.k_DR * D_const, 'o', 'color', 'white', 'markersize', mrksize - 3);
        end

        if strcmp(FP(cc, 2), 'Semistable')
            plot(FP(cc, 1), FP(cc, 1) + params_dict.k_DR * D_const, 'x', 'color', col, 'markersize', mrksize - 3);
        end

        if strcmp(FP(cc, 2), 'Stable')
            plot(FP(cc, 1), FP(cc, 1) + params_dict.k_DR * D_const, 'o', 'color', col, 'markersize', mrksize);
        end
    end

    % Plot threshold line
    threshold_position = params_dict.Theta / params_dict.k_BI;
    plot([threshold_position, threshold_position], [R_min, R_max], '--', 'color', 'red', 'linewidth', 3.0);

    % Set axis limits and labels
    xlim([B_min, B_max]);
    ylim([R_min, R_max]);
    xlabel('Extent of blood-brain barrier disruption B \approx I');
    ylabel('Degree of circuit remodeling R');
    title(['D = ' num2str(round(D_const * 100) / 100)]);

    % Save the figure
    %saveas(gcf, ['Figures/Fig7/Fig_7b' num2str(ii) '_inset.pdf']);

    % Close the figure
    %close(gcf);
end









%FUNCTION DEFINITION

%fun_dbdt function
function result_dbdt = fun_dbdt(x, params_dict, D_const)
    result_dbdt = -x + params_dict.k_IB * params_dict.k_BI * x + params_dict.K_SB * ...
        (exp(params_dict.k_IS * (params_dict.k_BI * x).^2 + ...
        params_dict.k_RS * params_dict.k_BR * x + params_dict.k_RS * ...
        params_dict.k_DR * D_const) - 1) / ...
        (exp(params_dict.k_IS * (params_dict.k_BI * x).^2 + ...
        params_dict.k_RS * params_dict.k_BR * x + params_dict.k_RS * ...
        params_dict.k_DR * D_const) + 1);
end

% fixedpointsfinder function 
function output = fixedpointsfinder(result_dbdt, a, b, n, params_dict, D_const)
    output = cell(0, 2);
    x = linspace(a, b, n+1);

    for ii = 2:(length(x)-1)
        if result_dbdt(x(ii), params_dict, D_const) <= 0
            if result_dbdt(x(ii-1), params_dict, D_const) > 0 && result_dbdt(x(ii+1), params_dict, D_const) < 0
                output{end+1, 1} = x(ii);
                output{end, 2} = 'Stable';
            elseif result_dbdt(x(ii-1), params_dict, D_const) < 0 && result_dbdt(x(ii+1), params_dict, D_const) > 0
                output{end+1, 1} = x(ii);
                output{end, 2} = 'Unstable';
            end
        elseif result_dbdt(x(ii), params_dict, D_const) >= 0
            if result_dbdt(x(ii-1), params_dict, D_const) < 0 && result_dbdt(x(ii+1), params_dict, D_const) < 0
                output{end+1, 1} = x(ii);
                output{end, 2} = 'Semistable';
            elseif result_dbdt(x(ii-1), params_dict, D_const) > 0 && result_dbdt(x(ii+1), params_dict, D_const) > 0
                output{end+1, 1} = x(ii);
                output{end, 2} = 'Semistable';
            end
        end
    end
end