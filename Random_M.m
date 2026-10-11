clc;
clear;
close all;

rng(1);
load('.\channel\common_scenario_K4.mat');

%% ================= System parameters =================
K = 4;
N = 4;
L = 5;
lambda = 0.1;
area_len = 8 * lambda;
d_min = lambda / 2;
noise_pow = 1e-9;
Pmax = 1;
PL0 = -40;
alpha_pathloss = 2.8;
BSPos = [0; 0];

Ms = [2, 4, 6, 8, 10];
Random_sumRate = zeros(length(Ms), 1);
numRealizations = 10;

user_pos = r_terminal.';
tx_pos_rel = d_tx;

%% ================= Antenna-number sweep =================
for idx = 1:length(Ms)
    M = Ms(idx);
    rate_tmp = zeros(numRealizations, 1);

    for r = 1:numRealizations
        ant_pos = generate_random_MA_positions(M, area_len, d_min);

        [rate_tmp(r), ~, ~, ~] = run_SDR_sca_sum_rate_ISAC( ...
            ant_pos, user_pos, tx_pos_rel, K, N, M, lambda, PL0, alpha_pathloss, ...
            noise_pow, Pmax, random_param, BSPos, L);

        fprintf('Random MA M=%d | realization %d/%d | Sum-Rate=%.4f bps/Hz\n', ...
            M, r, numRealizations, rate_tmp(r));
    end

    Random_sumRate(idx) = mean(rate_tmp);
end

save(strcat('.\res_M\Random_SDR.mat'), 'Ms', 'Random_sumRate');
