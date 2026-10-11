clc;
clear;
close all;

rng(2);
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
PL0 = -40;
alpha_pathloss = 2.8;
BSPos = [0; 0];

P_list = [0.2 0.4 0.6 0.8 1 1.2];
Random_sumRate = zeros(length(P_list), 1);
numRealizations = 10;

user_pos = r_terminal.';
tx_pos_rel = d_tx;

%% ================= Power sweep =================
for idx = 1:length(P_list)
    Pmax = P_list(idx);
    rate_tmp = zeros(numRealizations, 1);

    for r = 1:numRealizations
        ant_pos = generate_random_MA_positions(M, area_len, d_min);

        [rate_tmp(r), ~, ~, ~] = run_SDR_sca_sum_rate_ISAC( ...
            ant_pos, user_pos, tx_pos_rel, K, N, M, lambda, PL0, alpha_pathloss, ...
            noise_pow, Pmax, random_param, BSPos, L);

        fprintf('Random MA Pmax=%.1f | realization %d/%d | Sum-Rate=%.4f bps/Hz\n', ...
            Pmax, r, numRealizations, rate_tmp(r));
    end

    Random_sumRate(idx) = mean(rate_tmp);
end

save(strcat('.\res_P\Random_SDR_2.mat'), 'P_list', 'Random_sumRate');
