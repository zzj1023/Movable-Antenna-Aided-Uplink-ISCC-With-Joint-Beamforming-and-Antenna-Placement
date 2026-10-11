clc;
clear;

plot_dir = fileparts(mfilename('fullpath'));
result_dir = fullfile(plot_dir, '..', 'res_Area');

proposed_candidates = { ...
    'proposed_SDR_comput.mat', ...
    'proposed_SDR__comput.mat', ...
    'proposed_SDR.mat'};
candidate_exists = cellfun(@(name) ...
    isfile(fullfile(result_dir, name)), proposed_candidates);
candidate_idx = find(candidate_exists, 1);
if isempty(candidate_idx)
    error('No Proposed Area result file was found in %s.', result_dir);
end

proposed_data = load(fullfile(result_dir, ...
    proposed_candidates{candidate_idx}), 'totalsum_proposed');
ga_data = load(fullfile(result_dir, 'GA_SDR.mat'), 'GA_sumRate');
fpa_data = load(fullfile(result_dir, 'FPA_SDR.mat'), 'FPA_sumRate');
as_data = load(fullfile(result_dir, 'AS_SDR.mat'), 'finalRates_AS');
rpa_data = load(fullfile(result_dir, 'Random_SDR.mat'), 'totalsum_Random');

totalsum_proposed = proposed_data.totalsum_proposed(:);
totalsum_GA = ga_data.GA_sumRate(:);
totalsum_FPA = fpa_data.FPA_sumRate(:);
totalsum_AS = as_data.finalRates_AS(:);
totalsum_Random = rpa_data.totalsum_Random(:);

x = (3:8).';
expected_length = numel(x);
series_lengths = [numel(totalsum_proposed), numel(totalsum_GA), ...
    numel(totalsum_FPA), numel(totalsum_AS), numel(totalsum_Random)];
if any(series_lengths ~= expected_length)
    error('Every Area result must contain %d values; found [%s].', ...
        expected_length, num2str(series_lengths));
end

figure('Position', [100, 100, 720, 480]);
hold on;
set(gcf, 'DefaultTextFontName', 'Times New Roman', ...
    'DefaultTextInterpreter', 'tex');

h_proposed = plot(x, totalsum_proposed, '-ok', ...
    'MarkerFaceColor', 'w', 'LineWidth', 1);
h_GA = plot(x, totalsum_GA, '-->g', ...
    'MarkerFaceColor', 'w', 'LineWidth', 1);
h_FPA = plot(x, totalsum_FPA, '--^r', ...
    'MarkerFaceColor', 'w', 'LineWidth', 1);
h_AS = plot(x, totalsum_AS, '--db', ...
    'MarkerFaceColor', 'w', 'LineWidth', 1);
h_Random = plot(x, totalsum_Random, '--sm', ...
    'MarkerFaceColor', 'w', 'LineWidth', 1);

legend([h_proposed, h_GA, h_FPA, h_AS, h_Random], ...
    'Proposed', 'GA', 'FPA', 'AS', 'RPA', ...
    'Location', 'eastoutside', 'FontSize', 12);
xlabel('The normalized size of receive region, A/\lambda', 'FontSize', 14);
ylabel('Sum rate (bps/Hz)', 'FontSize', 14);
xticks(x);
xlim([x(1), x(end)]);
grid on;
set(gca, 'GridLineStyle', '--', 'FontName', 'Times New Roman');
