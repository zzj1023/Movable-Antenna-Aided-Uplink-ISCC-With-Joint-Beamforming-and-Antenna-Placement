clc;
clear;
close all;

proposed_data = load('..\res_K\proposed_SDR.mat');
ga_data = load('..\res_K\GA_SDR.mat');
fpa_data = load('..\res_K\FPA_SDR.mat');
as_data = load('..\res_K\AS_SDR.mat');
random_data = load('..\res_K\Random_SDR.mat');

% Older result files may contain only the performance vector. Use the K
% values from any result file that includes them.
all_data = {proposed_data, ga_data, fpa_data, as_data, random_data};
Ks = [];
for data_idx = 1:numel(all_data)
    if isfield(all_data{data_idx}, 'Ks')
        current_Ks = all_data{data_idx}.Ks(:);
        if isempty(Ks)
            Ks = current_Ks;
        elseif ~isequal(Ks, current_Ks)
            error('The K values in the result files are inconsistent.');
        end
    end
end

if isempty(Ks)
    error('None of the result files contains the variable Ks.');
end

totalsum_proposed = proposed_data.totalsum_proposed(:);
totalsum_GA = ga_data.GA_sumRate(:);
totalsum_FPA = fpa_data.FPA_sumRate(:);
totalsum_AS = as_data.finalRates_AS(:);
totalsum_Random = random_data.Random_sumRate(:);

result_lengths = [numel(totalsum_proposed), numel(totalsum_GA), ...
    numel(totalsum_FPA), numel(totalsum_AS), numel(totalsum_Random)];
if any(result_lengths ~= numel(Ks))
    error('Every performance vector must have the same length as Ks.');
end

figure;
set(gcf, 'DefaultTextFontName', 'Times New Roman', ...
    'DefaultTextInterpreter', 'tex');

plot(Ks, totalsum_proposed, '-k', 'LineWidth', 1);
hold on;
h_proposed = plot(Ks, totalsum_proposed, '-ok', ...
    'MarkerFaceColor', 'w', 'LineWidth', 1);

plot(Ks, totalsum_GA, '--g', 'LineWidth', 1);
h_GA = plot(Ks, totalsum_GA, '>g', ...
    'MarkerFaceColor', 'w', 'LineWidth', 1);

plot(Ks, totalsum_FPA, '--r', 'LineWidth', 1);
h_FPA = plot(Ks, totalsum_FPA, '^r', ...
    'MarkerFaceColor', 'w', 'LineWidth', 1);

plot(Ks, totalsum_AS, '--b', 'LineWidth', 1);
h_AS = plot(Ks, totalsum_AS, 'db', ...
    'MarkerFaceColor', 'w', 'LineWidth', 1);

plot(Ks, totalsum_Random, '--m', 'LineWidth', 1);
h_Random = plot(Ks, totalsum_Random, 'sm', ...
    'MarkerFaceColor', 'w', 'LineWidth', 1);

h_legend = legend( ...
    [h_proposed, h_GA, h_FPA, h_AS, h_Random], ...
    'Proposed', 'GA', 'FPA', 'AS', 'RPA', ...
    'Location', 'southeast');
set(h_legend, 'FontSize', 12);

xlabel('Number of devices, K', 'FontSize', 14);
ylabel('Sum-rate (bps/Hz)', 'FontSize', 14);
xticks(Ks);
xlim([min(Ks), max(Ks)]);
grid on;
ax = gca;
ax.GridLineStyle = '--';
