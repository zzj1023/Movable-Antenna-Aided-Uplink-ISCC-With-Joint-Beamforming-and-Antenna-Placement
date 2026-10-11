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
totalsum_proposed = zeros(length(area_list), 1);
best_pos_prev = [];
best_val_prev = 0;

user_pos = r_terminal.';
if ~exist('d_tx', 'var')
    d_tx = (0:N-1).' * (lambda / 2);
end
tx_pos_rel = d_tx;

%% ================= Area sweep =================
for idx = 1:length(area_list)
    area_len = area_list(idx);
    r_MA = generate_MA_positions(M, lambda, area_len);

    pso.I = 30;
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

    pos = pso.x_min + rand(pso.P, pso.dim) .* (pso.x_max - pso.x_min);
    if isempty(best_pos_prev)
        pos(1, :) = reshape(r_MA, [1, 2 * M]);
    else
        % The previous area's best solution is feasible in the larger nested area.
        pos(1, :) = best_pos_prev;
    end
    vel = pso.v_min + rand(pso.P, pso.dim) .* (pso.v_max - pso.v_min);

    pbest = pos;
    pbestVal = zeros(pso.P, 1);
    gbest = pbest(1, :);
    gbestVal = best_val_prev;

    for p = 1:pso.P
        ant_pos = reshape(pos(p, :), [2, M]);
        [sum_rate, ~, ~, ~] = run_SDR_sca_sum_rate_ISAC( ...
            ant_pos, user_pos, tx_pos_rel, K, N, M, lambda, PL0, alpha_pathloss, ...
            noise_pow, Pmax, random_param, BSPos, L);

        if isinf(sum_rate) || isnan(sum_rate) || sum_rate < 0
            continue;
        end

        pbestVal(p) = sum_rate;
        if sum_rate > gbestVal
            gbestVal = sum_rate;
            gbest = pos(p, :);
        end
    end

    for i = 1:pso.I
        for p = 1:pso.P
        w_weight = pso.w_max - (pso.w_max-pso.w_min)*i/pso.I;
        r1 = rand; r2 = rand;
        vel(p,:) = w_weight*vel(p,:) + ...
                   c1*r1*(pbest(p,:) - pos(p,:)) + ...
                   c2*r2*(gbest - pos(p,:));
        pos(p,:) = pos(p,:) + vel(p,:);
        % 限制范围
        pos(p,:) = min(max(pos(p,:), pso.rangeLB), pso.rangeUB);
        ant_pos = reshape(pos(p,:), [2,M]);
        pos(p,:) = ant_pos(:).';

        % ant_pos = enforce_min_distance(ant_pos, d_min, area_len);

        % 写回粒子
        % pos(p,:) = reshape(ant_pos,1,2*M);

          [sum_rate, ~, ~, ~] = run_SDR_sca_sum_rate_ISAC(ant_pos, user_pos, tx_pos_rel, K, N, M, lambda, PL0, alpha_pathloss, noise_pow, Pmax, random_param, BSPos, L);

            if isinf(sum_rate) || isnan(sum_rate) || sum_rate < 0
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

        fprintf('Area=%.1f lambda | PSO iter=%d/%d | Best Sum-Rate=%.4f bps/Hz\n', ...
            area_factors(idx), i, pso.I, gbestVal);
    end

    totalsum_proposed(idx) = gbestVal;
    best_pos_prev = gbest;
    best_val_prev = gbestVal;
end

if ~exist('.\res_Area', 'dir')
    mkdir('.\res_Area');
end
save(strcat('.\res_Area\proposed_SDR_comput.mat'), ...
    'area_factors', 'area_list', 'totalsum_proposed');
