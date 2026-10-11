clc;
clear;
close all;

Ks = [2, 4, 6, 8, 10];
N = 4;
M = 4;
L = 5;
lambda = 0.1;
area_len = 8 * lambda;
noise_pow = 1e-9;
Pmax = 1;
PL0 = -40;
alpha_pathloss = 2.8;
B = 10e6;
BSPos = [0; 0];

gridN = 5;
[xg, yg] = meshgrid(linspace(-area_len / 2, area_len / 2, gridN));
ant_pos_candidate = [xg(:).'; yg(:).'];
M_all = size(ant_pos_candidate, 2);
R = min(M + 1, M_all);

finalRates_AS = zeros(numel(Ks), 1);
maxIter = 3;

for idx = 1:numel(Ks)
    K = Ks(idx);
    env_file = fullfile('.\scenario_saved', ...
        ['K', num2str(K)], 'fixed_env.mat');
    env = load(env_file);
    user_pos = env.user_pos;
    tx_pos_rel = env.d_tx;
    current_param = env.random_param;
    Dk = env.Dk;
    Ck = env.Ck;
    fe = env.fe_k;
    computing = default_computing_config(K);
    computing.Lk = Dk;
    computing.Ck = Ck;
    computing.fe = fe;

    [~, ~, scores] = antenna_selection_energy( ...
        ant_pos_candidate, M, user_pos, tx_pos_rel, lambda, PL0, ...
        alpha_pathloss, current_param, K, N, BSPos, L);
    [~, score_order] = sort(scores, 'descend');
    topR = score_order(1:R);
    ant_pos_topR = ant_pos_candidate(:, topR);

    [~, sel_pos, ~] = antenna_selection_greedy( ...
        ant_pos_topR, M, user_pos, tx_pos_rel, K, N, lambda, PL0, ...
        alpha_pathloss, noise_pow, Pmax, B, Dk, Ck, fe, ...
        current_param, BSPos, L);

    totalRate_AS = zeros(maxIter, 1);
    for iteration = 1:maxIter
        totalRate_AS(iteration) = run_SDR_sca_sum_rate_ISAC( ...
            sel_pos, user_pos, tx_pos_rel, K, N, M, lambda, PL0, ...
            alpha_pathloss, noise_pow, Pmax, current_param, BSPos, L, ...
            computing);

        fprintf(['K=%d | AS evaluation=%d/%d | ', ...
            'Sum-rate=%.4f bps/Hz\n'], ...
            K, iteration, maxIter, totalRate_AS(iteration));
    end

    finalRates_AS(idx) = mean(totalRate_AS);
end

if ~exist('.\res_K', 'dir')
    mkdir('.\res_K');
end
save('.\res_K\AS_SDR.mat', 'Ks', 'finalRates_AS');
