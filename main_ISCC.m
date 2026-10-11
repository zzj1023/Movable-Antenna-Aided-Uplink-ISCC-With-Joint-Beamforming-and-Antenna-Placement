clc; 
clear; 
close all;

% 推荐开启并行计算池以加速运行 (如果您的电脑有多核CPU)
% 如果尚未开启，取消下面这行的注释
% if isempty(gcp('nocreate')), parpool; end

rng(1);
load('.\channel\common_scenario_K4.mat');
% load('.\channel\common_scenario_K6.mat');

%% ================= 系统参数 =================
K = 4;              % 设备数
N = 4;              % 每个设备发射天线数
M = 4;              % 基站侧可移动接收天线数 (MA)
L = 5;              % 多径数
lambda = 0.1;       % 波长 (m)
area_len = 8 * lambda;    % 可移动天线范围 ([-area_len/2, area_len/2])
d_min = lambda/2;          % 天线最小安全物理距离 (防碰撞)
noise_pow = 1e-9;         % 噪声功率
Pmax = 1;                  % 每个设备最大发射功率
% B = 10e6;                  % 信道带宽 Hz
PL0 = -40;                 % 参考路径损耗
alpha_pathloss = 2.8;

BSPos = [0;0];             % 基站参考中心位置

%% ================= 位置初始化 =================
% 这里假设您有一个 generate_MA_positions 函数，若没有请自己写一个随机生成不碰撞坐标的逻辑
r_MA = generate_MA_positions(M, lambda, area_len); 
user_pos = r_terminal.';
tx_pos_rel = d_tx;

%% ================= PSO 参数 =================
pso.I = 50;         % 外层最大迭代次数
pso.P = 70;         % 粒子群规模 (建议先改小测试，比如 20)
pso.dim = 2*M;      % 维度：每个天线有 x, y 两个坐标
pso.w_max = 0.9;
pso.w_min = 0.4;
c1 = 1.5; c2 = 1.5;

pso.rangeLB = -area_len/2;
pso.rangeUB = area_len/2;
pso.x_max = pso.rangeUB * ones(1, pso.dim);
pso.x_min = pso.rangeLB * ones(1, pso.dim);
pso.v_max = 0.04 * area_len * ones(1, pso.dim);
pso.v_min = -pso.v_max;

% 粒子结构初始化
pos = pso.x_min + rand(pso.P, pso.dim) .* (pso.x_max - pso.x_min);
pos(1,:) = reshape(r_MA, [1, 2*M]); % 把初始猜测放入第一个粒子
vel = pso.v_min + (rand(pso.P, pso.dim)) .* (pso.v_max - pso.v_min);

pbest = pos;
% 【修复1】：我们要最大化速率，因此历史最优初始值设为极小值 0
pbestVal = zeros(pso.P, 1); 

gbest = pbest(1,:);

gbestVal = 0; % 全局最优初始值
[~,bestIdx] = max(pbestVal);
gbest = pbest(bestIdx,:);
gbestVal = pbestVal(bestIdx);

totalsum_proposed = zeros(pso.I, 1);

%% ================= PSO 优化循环 =================
disp('================= 开始联合优化 =================');
disp('外层: PSO 寻找最佳天线位置');
disp('内层: WMMSE+SCA 求解最优发射波束与最大 Sum-Rate');

for i = 1:pso.I
    w_weight = pso.w_max - (pso.w_max - pso.w_min) * i / pso.I;
    
    % 注意：如果仿真太慢，强烈建议将这里的 for 改为 parfor (并行for循环)
    for p = 1:pso.P
        % 1. 速度和位置更新
        r1 = rand; r2 = rand;
        vel(p,:) = w_weight*vel(p,:) + ...
                   c1*r1*(pbest(p,:) - pos(p,:)) + ...
                   c2*r2*(gbest - pos(p,:));
        vel(p,:) = min( max(vel(p,:), pso.v_min), pso.v_max);
        pos(p,:) = pos(p,:) + vel(p,:);
        
        % 2. 限制在基站二维平面范围内
        pos(p,:) = min(max(pos(p,:), pso.rangeLB), pso.rangeUB);
        ant_pos = reshape(pos(p,:), [2, M]);
        pos(p,:) = ant_pos(:).';

        
       % [sum_rate, ~, ~, ~] = run_wmmse_sca_sum_rate_ISAC(ant_pos, user_pos, tx_pos_rel, K, N, M, lambda, PL0, alpha_pathloss, noise_pow, Pmax, random_param, BSPos, L);
       [sum_rate, rates, SINRs, Vcell] = run_SDR_sca_sum_rate_ISAC( ...
     ant_pos, user_pos, tx_pos_rel, K, N, M, lambda, PL0, alpha_pathloss, ...
     noise_pow, Pmax, random_param, BSPos, L);
    
   

  
        % 5. 【修复2】：更新个体最优 (寻找最大值，用 >)
        if sum_rate > pbestVal(p)
            pbestVal(p) = sum_rate;
            pbest(p,:) = pos(p,:);
        end

        % 更新全局最优 (寻找最大值，用 >)
        if sum_rate > gbestVal
            gbestVal = sum_rate;
            gbest = pos(p,:);
        end
end

    totalsum_proposed(i) = gbestVal;
    fprintf('>>> PSO 外层迭代 %d/%d：当前最大总速率 = %.4f bps/Hz \n', i, pso.I, gbestVal);
end

% save(strcat('.\sum_conve\proposed.mat'),'totalsum_proposed');
% save(strcat('.\sum_conve\proposed_1.mat'),'totalsum_proposed');
save(strcat('.\sum_conve\proposed_comput_2.mat'),'totalsum_proposed');

%% ================= 绘图 =================
figure;
plot(1:pso.I, totalsum_proposed, 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('PSO 迭代次数 (Iteration)', 'FontSize', 12); 
ylabel('系统总传输速率 (bps/Hz)', 'FontSize', 12);
title('ISAC-MA 系统 Sum-Rate 联合优化收敛曲线', 'FontSize', 14);
grid on;