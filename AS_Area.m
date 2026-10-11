clc; 
clear; 
close all;
rng(1);
load('.\channel\common_scenario_K4.mat');
% load('.\channel\common_scenario_K6.mat');
%% ================== 系统参数 ==================
K = 4;
M = 4;
N = 4;
lambda = 0.1;
PL0 = -40;          
alpha_pathloss = 2.8;
L = 5;  
noise_pow = 1e-9;
Pmax = 1;
B = 10e6;
Dk = 2e6 * ones(K,1);
Ck = 400e6 * ones(K,1);
fe = 3e9 * ones(K,1);
BSPos = [0;0];

%% ================= 用户位置（固定） =================

user_pos = r_terminal.';
tx_pos_rel = (0:N-1).' * (lambda/2);

%% ================= 候选天线 (5x5 网格) =================
area_len_list = [3 4 5 6 7 8] * lambda;
area_max = max(area_len_list);
grid_step = lambda/2;

x_vec = -area_max/2 : grid_step : area_max/2;
y_vec = -area_max/2 : grid_step : area_max/2;

[xg, yg] = meshgrid(x_vec, y_vec);
ant_pos_candidate_all = [xg(:)'; yg(:)'];

M_all = size(ant_pos_candidate_all, 2);

finalRates_AS = zeros(length(area_len_list), 1);        
fprintf('总候选天线数（最大area）：%d\n', M_all);
maxIter = 3;
for idx = 1:length(area_len_list)

    area_len = area_len_list(idx);
    valid_idx = find( ...
        abs(ant_pos_candidate_all(1,:)) <= area_len/2 & ...
        abs(ant_pos_candidate_all(2,:)) <= area_len/2 );

    ant_pos_candidate = ant_pos_candidate_all(:, valid_idx);
fprintf('候选天线数：%d\n', M_all);
    m_active = 4;

    %% ===== Step 1: Energy 选择前 R 个 =====
    % R = min(2*m_active, M_all);
    R = min(m_active + 2 , M_all);

  [~, ~, scores] = antenna_selection_energy( ...
    ant_pos_candidate, m_active, user_pos, tx_pos_rel, ...
    lambda, PL0, alpha_pathloss, random_param, K, N,BSPos,L);

    [~, idx_sort] = sort(scores, 'descend');
    topR = idx_sort(1:R);
    ant_pos_topR = ant_pos_candidate(:, topR);

    %% ===== Step 2: Greedy 精选 =====
  [sel_idx_local, sel_pos, delay_hist] = antenna_selection_greedy( ...
        ant_pos_topR, m_active, user_pos, tx_pos_rel, ...
        K, N, lambda, PL0, alpha_pathloss, ...
        noise_pow, Pmax, B, Dk, Ck, fe, random_param,BSPos,L);

    sel_idx = topR(sel_idx_local);
    fprintf('最终选中天线索引：');
    disp(sel_idx);

    %% ================== Step 3: Monte Carlo 评估 ==================

    totalRate_AS = zeros(maxIter,1);

    for I = 1:maxIter
    [totalRate_AS(I), ~, ~, ~] = run_SDR_sca_sum_rate_ISAC( ...
    sel_pos, user_pos, tx_pos_rel, K, N, m_active, lambda, PL0, alpha_pathloss, ...
    noise_pow, Pmax,random_param, BSPos, L);
        % totaldelay_ASall(I) = totaldelay_AS;

        fprintf('Monte Carlo %d/%d: Rate=%.4f\n', I, maxIter, totalRate_AS(I));
    end
    finalRates_AS(idx) = mean(totalRate_AS);
end
save(strcat('.\res_Area\AS_SDR.mat'), 'area_factors', 'area_list', 'finalRates_AS');
