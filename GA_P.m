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
area_len = 8 * lambda;
d_min = lambda / 2;
noise_pow = 1e-9;
PL0 = -40;
alpha_pathloss = 2.8;
BSPos = [0; 0];
P_list = [0.2 0.4 0.6 0.8 1 1.2];

user_pos = r_terminal.';
tx_pos_rel = d_tx;
computing = default_computing_config(K);

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

GA_sumRate = zeros(numel(P_list), 1);
GA_bestPos = cell(numel(P_list), 1);

for idx = 1:numel(P_list)
    Pmax = P_list(idx);
    [GA_sumRate(idx), GA_bestPos{idx}] = optimize_ma_ga( ...
        user_pos, tx_pos_rel, K, N, M, lambda, PL0, alpha_pathloss, ...
        noise_pow, Pmax, random_param, BSPos, L, area_len, d_min, ...
        ga, computing, []);
    fprintf('Pmax=%.1f | GA sum-rate=%.4f bps/Hz\n', ...
        Pmax, GA_sumRate(idx));
end

if ~exist('.\res_P', 'dir')
    mkdir('.\res_P');
end
save('.\res_P\GA_SDR_2.mat', 'P_list', 'GA_sumRate', 'GA_bestPos', 'ga');
% save('.\res_P\GA_SDR_comput.mat', 'P_list', 'GA_sumRate', 'GA_bestPos', 'ga');