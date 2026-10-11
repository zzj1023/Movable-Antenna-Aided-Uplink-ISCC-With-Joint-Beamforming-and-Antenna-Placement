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
d_min = lambda / 2;
noise_pow = 1e-9;
Pmax = 1;
PL0 = -40;
alpha_pathloss = 2.8;
BSPos = [0; 0];

area_factors = [3 4 5 6 7 8];
area_list = area_factors * lambda;
Random_sumRate = zeros(length(area_list), 1);
numRealizations = 10;

user_pos = r_terminal.';
if ~exist('d_tx', 'var')
    d_tx = (0:N-1).' * (lambda / 2);
end
tx_pos_rel = d_tx;

%% ================= Area sweep =================
for idx = 1:length(area_list)
    area_len = area_list(idx);
    rate_tmp = zeros(numRealizations, 1);

    for r = 1:numRealizations
        ant_pos = generate_random_MA_positions(M, area_len, d_min);

        [rate_tmp(r), ~, ~, ~] = run_SDR_sca_sum_rate_ISAC( ...
            ant_pos, user_pos, tx_pos_rel, K, N, M, lambda, PL0, alpha_pathloss, ...
            noise_pow, Pmax, random_param, BSPos, L);

        fprintf('Random MA Area=%.1f lambda | realization %d/%d | Sum-Rate=%.4f bps/Hz\n', ...
            area_factors(idx), r, numRealizations, rate_tmp(r));
    end

    Random_sumRate(idx) = mean(rate_tmp);
end

if ~exist('.\res_Area', 'dir')
    mkdir('.\res_Area');
end
save(strcat('.\res_Area\Random_SDR.mat'), 'area_factors', 'area_list', 'Random_sumRate');
