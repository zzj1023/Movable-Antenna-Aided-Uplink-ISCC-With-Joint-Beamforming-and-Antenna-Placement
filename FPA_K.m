clc;
clear;
close all;
rng(1);

Ks = [2, 4, 6, 8, 10];
N = 4;
M = 4;
L = 5;
lambda = 0.1;
noise_pow = 1e-9;
Pmax = 1;
PL0 = -40;
alpha_pathloss = 2.8;
BSPos = [0; 0];

FPA_sumRate = zeros(numel(Ks), 1);
maxIter = 1;

for idx = 1:numel(Ks)
    K = Ks(idx);
    env_file = fullfile('.\scenario_saved', ...
        ['K', num2str(K)], 'fixed_env.mat');
    env = load(env_file);

    user_pos = env.user_pos;
    tx_pos_rel = env.d_tx;
    random_param = env.random_param;
    ant_pos = env.ant_pos;
    computing = default_computing_config(K);
    computing.Lk = env.Dk;
    computing.Ck = env.Ck;
    computing.fe = env.fe_k;

    totalSR_FPA = zeros(maxIter, 1);
    for iteration = 1:maxIter
        totalSR_FPA(iteration) = run_SDR_sca_sum_rate_ISAC( ...
            ant_pos, user_pos, tx_pos_rel, K, N, M, lambda, PL0, ...
            alpha_pathloss, noise_pow, Pmax, random_param, BSPos, L, ...
            computing);
    end
    FPA_sumRate(idx) = min(totalSR_FPA);

    fprintf('K=%d | FPA sum-rate=%.4f bps/Hz\n', ...
        K, FPA_sumRate(idx));
end

if ~exist('.\res_K', 'dir')
    mkdir('.\res_K');
end
save('.\res_K\FPA_SDR.mat', 'Ks', 'FPA_sumRate');
