clc;
clear;
close all;
rng(1);

%% System parameters
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
BSPos = [0; 0];

totalsum_proposed = zeros(numel(Ks), 1);

for idx = 1:numel(Ks)
    K = Ks(idx);
    env_file = fullfile('.\scenario_saved', ...
        ['K', num2str(K)], 'fixed_env.mat');
    env = load(env_file);
    user_pos = env.user_pos;
    tx_pos_rel = env.d_tx;
    random_param = env.random_param;
    computing = default_computing_config(K);
    computing.Lk = env.Dk;
    computing.Ck = env.Ck;
    computing.fe = env.fe_k;

    % Same MA initialization used by the M/P experiments.
    r_MA = generate_MA_positions(M, lambda, area_len);

    %% PSO parameters
    pso.I = 50;
    pso.P = 70;
    pso.dim = 2 * M;
    pso.w_max = 0.9;
    pso.w_min = 0.4;
    c1 = 1.5;
    c2 = 1.5;
    pso.rangeLB = -area_len / 2;
    pso.rangeUB = area_len / 2;
    pso.x_max = pso.rangeUB * ones(1, pso.dim);
    pso.x_min = pso.rangeLB * ones(1, pso.dim);
    pso.v_max = 0.04 * area_len * ones(1, pso.dim);
    pso.v_min = -pso.v_max;

    pos = pso.x_min + rand(pso.P, pso.dim) .* ...
        (pso.x_max - pso.x_min);
    pos(1, :) = reshape(r_MA, [1, 2 * M]);
    vel = pso.v_min + rand(pso.P, pso.dim) .* ...
        (pso.v_max - pso.v_min);

    pbest = pos;
    pbestVal = zeros(pso.P, 1);
    gbest = pbest(1, :);
    gbestVal = 0;

    %% PSO optimization
    for iteration = 1:pso.I
        w_weight = pso.w_max - ...
            (pso.w_max - pso.w_min) * iteration / pso.I;

        for p = 1:pso.P
            r1 = rand;
            r2 = rand;
            vel(p, :) = w_weight * vel(p, :) + ...
                c1 * r1 * (pbest(p, :) - pos(p, :)) + ...
                c2 * r2 * (gbest - pos(p, :));
            pos(p, :) = pos(p, :) + vel(p, :);
            pos(p, :) = min(max(pos(p, :), ...
                pso.rangeLB), pso.rangeUB);

            ant_pos = reshape(pos(p, :), [2, M]);
            pos(p, :) = ant_pos(:).';

            sum_rate = run_SDR_sca_sum_rate_ISAC( ...
                ant_pos, user_pos, tx_pos_rel, K, N, M, lambda, ...
                PL0, alpha_pathloss, noise_pow, Pmax, random_param, ...
                BSPos, L, computing);

            if ~isfinite(sum_rate) || sum_rate < 0
                continue;
            end

            if sum_rate > pbestVal(p)
                pbestVal(p) = sum_rate;
                pbest(p, :) = pos(p, :);
            end
            if sum_rate > gbestVal
                gbestVal = sum_rate;
                gbest = pos(p, :);
            end
        end

        fprintf(['K=%d | PSO iteration=%d/%d | ', ...
            'Best sum-rate=%.4f bps/Hz\n'], ...
            K, iteration, pso.I, gbestVal);
    end

    totalsum_proposed(idx) = gbestVal;
end

if ~exist('.\res_K', 'dir')
    mkdir('.\res_K');
end
save('.\res_K\proposed_SDR.mat', ...
    'Ks', 'totalsum_proposed', 'pso');
