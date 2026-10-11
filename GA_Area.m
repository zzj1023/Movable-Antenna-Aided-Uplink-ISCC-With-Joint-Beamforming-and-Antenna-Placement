clc;
clear;
close all;
rng(1);
load('.\channel\common_scenario_K4.mat');

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

user_pos = r_terminal.';
tx_pos_rel = d_tx;
computing = default_computing_config(K);
GA_sumRate = zeros(numel(area_list), 1);
GA_bestPos = cell(numel(area_list), 1);
initial_pos = [];

for idx = 1:numel(area_list)
    area_len = area_list(idx);
    ga.numGenerations = 30;
    ga.populationSize = 70;
    ga.dim = 2 * M;
    ga.crossoverRate = 0.9;
    ga.mutationRate = 1 / ga.dim;
    ga.eliteCount = 2;
    ga.tournamentSize = 3;
    ga.blendAlpha = 0.5;
    ga.sigmaMax = 0.10 * area_len;
    ga.sigmaMin = 0.01 * area_len;

    [GA_sumRate(idx), GA_bestPos{idx}] = optimize_ma_ga( ...
        user_pos, tx_pos_rel, K, N, M, lambda, PL0, alpha_pathloss, ...
        noise_pow, Pmax, random_param, BSPos, L, area_len, d_min, ...
        ga, computing, initial_pos);
    initial_pos = GA_bestPos{idx};
    fprintf('Area=%.1f lambda | GA sum-rate=%.4f bps/Hz\n', ...
        area_factors(idx), GA_sumRate(idx));
end

if ~exist('.\res_Area', 'dir')
    mkdir('.\res_Area');
end
save('.\res_Area\GA_SDR.mat', ...
    'area_factors', 'area_list', 'GA_sumRate', 'GA_bestPos');
