clc;
clear;
close all;

rng(1);
load('.\channel\common_scenario_K4.mat');

%% ================= System parameters =================
K = 4;
N = 4;
M = 4;
L = 5;
lambda = 0.1;
area_len = 8 * lambda;
d_min = lambda / 2;
noise_pow = 1e-9;
Pmax = 1;
PL0 = -40;
alpha_pathloss = 2.8;
BSPos = [0; 0];

user_pos = r_terminal.';
tx_pos_rel = d_tx;

%% ================= Random MA baseline =================
maxIter = 50;
random_samples = zeros(maxIter, 1);

for I = 1:maxIter
    ant_pos = generate_random_MA_positions(M, area_len, d_min);

    [sum_rate, ~, ~, ~] = run_SDR_sca_sum_rate_ISAC( ...
        ant_pos, user_pos, tx_pos_rel, K, N, M, lambda, PL0, alpha_pathloss, ...
        noise_pow, Pmax, random_param, BSPos, L);

    random_samples(I) = sum_rate;
    fprintf('Random MA Iter %d/%d: Sum-Rate = %.4f bps/Hz\n', I, maxIter, sum_rate);
end

random_running_mean = cumsum(random_samples) ./ (1:maxIter).';
random_mean_rate = mean(random_samples);

% RPA has no iterative position optimization. Use its Monte Carlo mean as
% the constant reference in the convergence comparison.
totalsum_Random = random_mean_rate * ones(maxIter, 1);
save(strcat('.\sum_conve\Random_comput.mat'), 'totalsum_Random', ...
    'random_samples', 'random_running_mean', 'random_mean_rate');

figure;
plot(1:maxIter, random_running_mean, '-xm', ...
    'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Random realization index', 'FontSize', 12);
ylabel('Sum-Rate (bps/Hz)', 'FontSize', 12);
title('RPA Monte Carlo Mean', 'FontSize', 14);
grid on;
