clc; 
clear; 
close all;
load('.\channel\common_scenario_K4.mat');
% load('.\channel\common_scenario_K6.mat');
%% ================== 系统参数 ==================
K = 4;    
N = 4;
L = 5;
lambda = 0.1;
PL0 = -40;          
alpha_pathloss = 2.8;

noise_pow = 1e-9;
Pmax = 1;
B = 10e6;

Dk = 2e6 * ones(K,1);
Ck = 400e6 * ones(K,1);
fe = 3e9 * ones(K,1);
BSPos = [0;0];


%% ================= 用户位置（固定） =================
% r_terminal = [30*rand(1,K); 30*rand(1,K)];
% r_MA = generate_MA_positions(M, lambda, area_len);
user_pos = r_terminal.';
tx_pos_rel = (0:N-1).' * (lambda/2);

%% ================= 候选天线 (5x5 网格) =================
area_len = 8 * lambda;
gridN = 5;

[xg, yg] = meshgrid(linspace(-area_len*0.5, area_len*0.5, gridN));
ant_pos_candidate = [xg(:)'; yg(:)'];
M_all = size(ant_pos_candidate, 2);

fprintf('候选天线数：%d\n', M_all);

%% ================= 天线选择目标数 =================
m_active = 4;

%% ===== Step 1: Energy 选择前 R 个 =====
R = min(2*m_active, M_all);
% R = min(m_active + 2, M_all);

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
maxIter = 50;
totalSR_AS = zeros(maxIter,1);

for I = 1:maxIter
    % [totalSR_AS(I), ~, ~, ~] = run_wmmse_sca_sum_rate_ISAC( ...
    % sel_pos, user_pos, tx_pos_rel, K, N, m_active, lambda, PL0, alpha_pathloss, ...
    % noise_pow, Pmax,random_param, BSPos, L);
     [totalSR_AS(I), ~, ~, ~] = run_SDR_sca_sum_rate_ISAC( ...
     sel_pos, user_pos, tx_pos_rel, K, N, m_active, lambda, PL0, alpha_pathloss, ...
     noise_pow, Pmax, random_param, BSPos, L);
    


    fprintf('Monte Carlo %d/%d: 总速率 = %.4f bps/Hz \n', I, maxIter, totalSR_AS(I));
end

% 保存结果
%% ================================
save('.\sum_conve\AS_comput.mat', 'totalSR_AS');
% save('.\Iter_50\AS_selection.mat', 'totalDelay_AS');
% save('.\channel\AS_init_pos.mat', 'sel_pos');
%% ================== 绘图 ==================

