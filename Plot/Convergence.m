clc;
clear;

plot_dir = fileparts(mfilename('fullpath'));
result_dir = fullfile(plot_dir, '..', 'sum_conve');

proposed_data = load(fullfile(result_dir, 'proposed_comput.mat'), ...
    'totalsum_proposed');
ga_data = load(fullfile(result_dir, 'GA_1.mat'), 'totalsum_GA');
fpa_data = load(fullfile(result_dir, 'FPA_comput.mat'), 'totalsum_FPA');
as_data = load(fullfile(result_dir, 'AS_comput.mat'), 'totalSR_AS');
random_data = load(fullfile(result_dir, 'Random_comput.mat'));

totalsum_proposed = proposed_data.totalsum_proposed(:);
totalsum_GA = ga_data.totalsum_GA(:);
totalsum_FPA = fpa_data.totalsum_FPA(:);
totalsum_AS = as_data.totalSR_AS(:);

if isfield(random_data, 'random_mean_rate')
    random_mean_rate = random_data.random_mean_rate;
elseif isfield(random_data, 'random_samples')
    random_mean_rate = mean(random_data.random_samples(:));
elseif isfield(random_data, 'totalsum_Random')
    random_mean_rate = mean(random_data.totalsum_Random(:));
else
    error('Random_comput.mat does not contain recognized RPA results.');
end

n = min([50, numel(totalsum_proposed), numel(totalsum_GA), ...
    numel(totalsum_FPA), numel(totalsum_AS)]);
if n < 1
    error('The convergence result files contain no samples.');
end

x = (1:n).';
totalsum_Random = random_mean_rate * ones(n, 1);
marker_step = 5;
marker_idx = unique([1:marker_step:n, n]);

figure;
hold on;
set(gcf, 'DefaultTextFontName', 'Times New Roman', ...
    'DefaultTextInterpreter', 'tex');

h_proposed = plot(x, totalsum_proposed(1:n), '-ok', ...
    'MarkerIndices', marker_idx, 'MarkerFaceColor', 'w', 'LineWidth', 1);
h_GA = plot(x, totalsum_GA(1:n), '-->g', ...
    'MarkerIndices', marker_idx, 'MarkerFaceColor', 'w', 'LineWidth', 1);
h_FPA = plot(x, totalsum_FPA(1:n), '--^r', ...
    'MarkerIndices', marker_idx, 'MarkerFaceColor', 'w', 'LineWidth', 1);
h_AS = plot(x, totalsum_AS(1:n), '--db', ...
    'MarkerIndices', marker_idx, 'MarkerFaceColor', 'w', 'LineWidth', 1);
h_Random = plot(x, totalsum_Random, '--sm', ...
    'MarkerIndices', marker_idx, 'MarkerFaceColor', 'w', 'LineWidth', 1);

legend([h_proposed, h_GA, h_FPA, h_AS, h_Random], ...
    'Proposed', 'GA', 'FPA', 'AS', 'RPA', ...
    'Location', 'southeast', 'FontSize', 12);
xlabel('Number of iterations', 'FontSize', 14);
ylabel('Sum rate (bps/Hz)', 'FontSize', 14);
xticks(0:5:n);
xlim([1, n]);
grid on;
set(gca, 'GridLineStyle', '--', 'FontName', 'Times New Roman');
