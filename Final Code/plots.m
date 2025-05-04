clear all; close all; clc;
addpath('data')

%% Plot 1: RMSE Comparison from 64x64 to 256x256 Spectral
% Load RMSE data for all methods
load('fvm_ic1_64x64.mat'); 
load('fdm_ic1_64x64.mat'); 
load('spectral_ic1_64x64.mat');

% Time points
save_times = 0:4:32;

% Plot
figure;
plot(save_times, rmse_fvm, 's-', 'LineWidth', 1.5, 'DisplayName', 'FVM');
hold on;
plot(save_times, rmse_fdm, 'o-', 'LineWidth', 1.5, 'DisplayName', 'FDM');
plot(save_times, rmse_spectral, '^-', 'LineWidth', 1.5, 'DisplayName', 'Spectral');

xlabel('Time (s)', 'FontSize', 12);
ylabel('RMSE', 'FontSize', 12);
title('RMSE Over Time for FVM, FDM, and Spectral Methods', 'FontSize', 14);
legend('Location', 'northwest', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 12);

%% Plot 2: Error Growth Rate (dRMSE/dt)
% Compute error growth rates using finite differences
growth_rate_fvm = diff(rmse_fvm)/4;        % Time step Δt = 4 seconds
growth_rate_fdm = diff(rmse_fdm)/4;
growth_rate_spectral = diff(rmse_spectral)/4;
growth_times = save_times(2:end);           % Times for growth rate data

% Plot settings
figure;
plot(growth_times, growth_rate_fvm, 's-', 'LineWidth', 1.5, 'DisplayName', 'FVM');
hold on;
plot(growth_times, growth_rate_fdm, 'o-', 'LineWidth', 1.5, 'DisplayName', 'FDM');
plot(growth_times, growth_rate_spectral, '^-', 'LineWidth', 1.5, 'DisplayName', 'Spectral');

xlabel('Time (s)', 'FontSize', 12);
ylabel('Error Growth Rate (RMSE/s)', 'FontSize', 12);
title('Error Growth Rate Comparison', 'FontSize', 14);
legend('Location', 'northwest', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 12);

%% Plot 3: Timing Comparison with Horizontal Value Labels

% Load timing data
load('fvm_timing_results.mat'); 
load('fdm_timing_results.mat'); 
load('spectral_timing_results.mat');

% Organize data
methods = {'FVM', 'FDM', 'Spectral'};
initial_conditions = {'Single Vortex', 'Colliding Vortices', 'Multiple Vortices', 'Dipole'};
timing_data = [fvm_timing_results; fdm_timing_results; spectral_timing_results];

% Create figure
figure;
h = bar(timing_data', 'grouped');
set(gca, 'XTickLabel', initial_conditions, 'FontSize', 12, 'YTick', 0:5:max(timing_data(:))+10);
ylabel('Computational Time (s)', 'FontSize', 12);
title('Method Comparison: Time per Initial Condition', 'FontSize', 14);
legend(methods, 'Location', 'northwest', 'FontSize', 10);
grid on;

% Customize colors
colors = [0.2 0.6 0.8;   % FVM (blue)
          0.8 0.3 0.3;   % FDM (red)
          0.4 0.8 0.4];  % Spectral (green)
for m = 1:3
    h(m).FaceColor = colors(m,:);
end

% Calculate bar positions
bar_width = 0.8;
group_width = bar_width/3;
x_positions = zeros(4,3);  % [IC × Method]

for ic = 1:4
    base_x = ic - bar_width/2;
    for meth = 1:3
        x_positions(ic,meth) = base_x + (meth-1)*group_width + group_width/2;
    end
end

% Add horizontal text labels
for ic = 1:4
    for meth = 1:3
        bar_height = timing_data(meth,ic);
        text(x_positions(ic,meth), bar_height + 0.5,...
            sprintf('%.1f', bar_height),...
            'Rotation', 0,...
            'FontSize', 10,...
            'HorizontalAlignment', 'center',...
            'VerticalAlignment', 'bottom');
    end
end



%% Plot 4: Final Error Comparison (Bar Chart)
figure;
methods = {'FVM', 'FDM', 'Spectral'};
final_errors = [rmse_fvm(end), rmse_fdm(end), rmse_spectral(end)];
bar(final_errors);
set(gca, 'XTickLabel', methods, 'FontSize', 12);
ylabel('Final RMSE at t=32s', 'FontSize', 12);
title('Final Error Comparison', 'FontSize', 14);
grid on;


%% Plot 5: Time against Resolution (Improved Readability)

resolutions = [32, 64, 128];
methods = {'FVM', 'FDM', 'Spectral'};

% Preallocate time matrix: rows = resolutions, columns = methods
time_matrix = zeros(length(resolutions), length(methods));

% Load data into matrix
for meth_idx = 1:length(methods)
    method = methods{meth_idx};
    for res_idx = 1:length(resolutions)
        m = resolutions(res_idx);
        time_file = sprintf('time_%s_%d.mat', lower(method), m);
        data = load(time_file);
        vars = fieldnames(data);
        time_matrix(res_idx, meth_idx) = data.(vars{1});
    end
end

% Create figure with grouped bars
figure;
hb = bar(time_matrix, 'grouped');
set(gca, 'YScale', 'log', 'FontSize', 12);
grid on;

% Axis labels and title
xlabel('Grid Resolution', 'FontSize', 12);
ylabel('Computational Time (s)', 'FontSize', 12);
title('Computational Time vs Resolution', 'FontSize', 14);
xticklabels({'32²', '64²', '128²'});
legend(methods, 'Location', 'northwest', 'FontSize', 10);

% Add centered time labels with black text
for i = 1:length(methods)
    hb(i).FaceColor = colors(i,:);
    x_positions = hb(i).XEndPoints; % Get actual bar centers
    
    for j = 1:length(resolutions)
        text(x_positions(j), time_matrix(j,i)*1.05, ... % 5% above bar
            sprintf('%.1f', time_matrix(j,i)), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', 10, ...
            'Color', 'k'); % Changed to black text
    end
end

% Set y-axis limits with 20% headroom
max_time = max(time_matrix(:));
ylim([0.1, max_time*1.5]); % Logarithmic y-axis