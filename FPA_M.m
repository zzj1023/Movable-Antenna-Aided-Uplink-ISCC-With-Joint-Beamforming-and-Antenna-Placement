clc; 
clear; 
close all;
rng(1);
load('.\channel\common_scenario_K4.mat');
%% ================= 系统参数 =================
K = 4;              % 设备数
N = 4;              % 每个设备天线数
% M = 4;              % 基站天线数（固定）
L = 5;              % 多径数
lambda = 0.1;       % 波长 (m)
noise_pow = 1e-9;   % 噪声功率
Pmax = 1;           % 每个设备最大发射功率
% B = 10e6;            % 信道带宽 Hz
PL0 = -40;            % 参考路径损耗
alpha_pathloss = 2.8;

Dk = 2e6 * ones(K,1);   % 上传数据量 (bit)
Ck = 400e6 * ones(K,1); % 任务计算量 (CPU cycles)
fe = 3e9 * ones(K,1);   % MEC CPU分配 Hz
BSPos = [0;0];

%% ================= 位置初始化 =================

% 设备位置 (随机分布)
% r_terminal = [30*rand(1,K); 30*rand(1,K)];
% 设备天线线性排列坐标
d_tx = (0:N-1).' * (lambda/2);
Ms = [2,4,6,8,10];
FPA_sumRate = zeros(length(Ms),1); 
%% ================= 迭代次数设置 =================
maxIter = 1;  % 可以模拟多次随机初始化收敛效果
for idx = 1:length(Ms)
    M = Ms(idx);

    area_len = 8 * lambda;
    d = lambda/2;                % 理想间隔（0.05 m）
    n_side = ceil(sqrt(M));      % 阵列每边天线数

    [xg, yg] = meshgrid(0:n_side-1, 0:n_side-1);
    xg = xg(:)*d;
    yg = yg(:)*d;

    % 取前M个天线
    ant_pos = [xg(1:M), yg(1:M)].';    % 2×M

    % 阵列居中到原点
    ant_pos = ant_pos - mean(ant_pos,2);

    totalSR_FPA = zeros(maxIter,1);

%% ================= 主循环: 固定天线布局 =================
for I = 1:maxIter
    user_pos = r_terminal.';   % 设备位置
    tx_pos_rel = d_tx;
     [H, A, H_int] = channel_model_ISAC_deterministic( ...
    user_pos, ant_pos, tx_pos_rel, lambda, PL0, alpha_pathloss, ...
    random_param, K, BSPos, L);
      V_init = cell(K,1);
            for k = 1:K
                heq = sum(H(:,:,k), 1).'; % 简单合并信道
                if norm(heq) > 0
                    V_init{k} = heq / norm(heq) * sqrt(Pmax);
                else
                    tmp = randn(N,1)+1j*randn(N,1);
                    V_init{k} = tmp/norm(tmp)*sqrt(Pmax);
                end
            end
      [sum_rate, rates, SINRs, Vcell] = run_SDR_sca_sum_rate_ISAC( ...
    ant_pos, user_pos, tx_pos_rel, K, N, M, lambda, PL0, alpha_pathloss, ...
    noise_pow, Pmax,random_param, BSPos, L);

    totalSR_FPA(I) = sum_rate;

    fprintf('M=%d 固定天线 第 %d/%d 次：总Rate = %.4f bps\n', ...
            M, I, maxIter, totalSR_FPA);
    % fprintf('固定天线方案迭代 %d/%d：总时延 = %.4f 秒\n', I, maxIter, totalSR_FPA);
end
    FPA_sumRate(idx) = min(totalSR_FPA);

    % 保存每次迭代曲线（可选）
    % Delay_record{idx} = totalDelay_FPA;
end
%% ================= 绘图 =================
figure;
plot(Ms, FPA_sumRate, 'r-s', 'LineWidth', 1.8);
xlabel('基站天线数 M'); ylabel('最优总时延 (s)');
title('固定天线布局：不同 M 时的最优时延');
grid on;

%% ================= 保存结果 =================
% save(strcat('.\res_M\FPA.mat'), 'FPA_sumRate');
save(strcat('.\res_M\FPA_SDR.mat'), 'FPA_sumRate');