clc;
clear;
close all;
rng(1);
load('.\channel\common_scenario_K4.mat');

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
Ms = [2 4 6 8 10];

user_pos = r_terminal.';
tx_pos_rel = d_tx;
computing = default_computing_config(K);
GA_sumRate = zeros(numel(Ms), 1);
GA_bestPos = cell(numel(Ms), 1);

for idx = 1:numel(Ms)
    M = Ms(idx);
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
        ga, computing, []);
    fprintf('M=%d | GA sum-rate=%.4f bps/Hz\n', M, GA_sumRate(idx));
end

if ~exist('.\res_M', 'dir')
    mkdir('.\res_M');
end
save('.\res_M\GA_SDR.mat', 'Ms', 'GA_sumRate', 'GA_bestPos');
