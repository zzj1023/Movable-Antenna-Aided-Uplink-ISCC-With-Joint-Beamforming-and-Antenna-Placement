clc;
clear;
close all;
rng(1);

Ks = [2, 4, 6, 8, 10];
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
numRealizations = 15;

Random_sumRate = zeros(numel(Ks), 1);

for idx = 1:numel(Ks)
    K = Ks(idx);
    env_file = fullfile('.\scenario_saved', ...
        ['K', num2str(K)], 'fixed_env.mat');
    env = load(env_file);
    user_pos = env.user_pos;
    tx_pos_rel = env.d_tx;
    current_param = env.random_param;
    computing = default_computing_config(K);
    computing.Lk = env.Dk;
    computing.Ck = env.Ck;
    computing.fe = env.fe_k;
    rate_tmp = zeros(numRealizations, 1);

    for realization = 1:numRealizations
        ant_pos = generate_random_MA_positions(M, area_len, d_min);
        rate_tmp(realization) = run_SDR_sca_sum_rate_ISAC( ...
            ant_pos, user_pos, tx_pos_rel, K, N, M, lambda, PL0, ...
            alpha_pathloss, noise_pow, Pmax, current_param, BSPos, L, ...
            computing);

        fprintf(['K=%d | RPA realization=%d/%d | ', ...
            'Sum-rate=%.4f bps/Hz\n'], K, realization, ...
            numRealizations, rate_tmp(realization));
    end

    Random_sumRate(idx) = mean(rate_tmp);
end

if ~exist('.\res_K', 'dir')
    mkdir('.\res_K');
end
save('.\res_K\Random_SDR.mat', 'Ks', 'Random_sumRate');
