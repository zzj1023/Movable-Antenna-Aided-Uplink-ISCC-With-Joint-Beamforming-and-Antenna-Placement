clc; 
clear; 
close all;
rng(1);
load('.\channel\common_scenario_K4.mat');
% load('.\channel\common_scenario_K6.mat');
%% ================= 系统参数 =================
K = 4;              % 设备数
N = 4;              % 每个设备天线数
M = 4;              % 可移动接收天线数
L = 5;              % 多径数
lambda = 0.1;       % 波长 (m)
area_len = 8 * lambda;    % 可移动天线范围
d_min = lambda/2;
noise_pow = 1e-9;   % 噪声功率
% Pmax = 1;         % 每个设备最大发射功率
% B = 10e6;            % 信道带宽 Hz
PL0 = -40;            % 参考路径损耗
alpha_pathloss = 2.8;
BSPos = [0;0];

%% ================= 位置初始化 =================
P_list = [0.2 0.4 0.6 0.8 1 1.2];                         % 不同天线数
totalsum_proposed = zeros(length(P_list),1);   % 每个 M 的最终最优时延

%% ================= 主循环 =================
for idx = 1:length(P_list)
    Pmax = P_list(idx);

    % 基站天线初始坐标 (二维)
    r_MA = generate_MA_positions(M, lambda, area_len);

    % 设备天线线性排列坐标
    d_tx = (0:N-1).' * (lambda/2);
    user_pos = r_terminal.';
    tx_pos_rel = d_tx;

    %% ================= PSO 参数 =================
    pso.I = 30;
    pso.P = 70;
    pso.dim = 2*M;
    pso.w_max = 0.9;
    pso.w_min = 0.4;
    c1 = 1.5; c2 = 1.5;
    pso.rangeLB = -area_len/2;
    pso.rangeUB = area_len/2;
    pso.x_max = pso.rangeUB*ones(1,pso.dim);
    pso.x_min = pso.rangeLB*ones(1,pso.dim);
    pso.v_max = 0.04*area_len*ones(1,pso.dim);
    pso.v_min = -pso.v_max;
    % 粒子结构初始化
    pos = pso.x_min + rand(pso.P,pso.dim).*(pso.x_max-pso.x_min);
    pos(1,:) = reshape(r_MA, [1, 2*M]);
    vel = pso.v_min + (rand(pso.P,pso.dim)).*(pso.v_max-pso.v_min);

    % pos = area_len * rand(numParticles, 2*M);
    % vel = zeros(numParticles, 2*M);
    pbest = pos;
    % 【修复1】：我们要最大化速率，因此历史最优初始值设为极小值 0
    pbestVal = zeros(pso.P, 1); 

    gbest = pbest(1,:);

    gbestVal = 0; % 全局最优初始值

    [~,bestIdx] = max(pbestVal);
    gbest = pbest(bestIdx,:);
    gbestVal = pbestVal(bestIdx);

    % totalDelay_proposed = zeros(5,1);

%% ================= PSO 优化循环 =================
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
  
        % 过滤异常值
        if isinf(sum_rate) || isnan(sum_rate) || sum_rate < 0
            fprintf('警告：粒子 %d 的sum_rate异常（%.4f），跳过更新\n', p, sum_rate);
            continue;
        end

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

    % totalDelay_proposed(i) = gbestVal;
           fprintf('Pmax=%d | PSO迭代=%d/%d | 最优和速率=%.4f bps\n', ...
            Pmax, i, pso.I, gbestVal);
    end
    totalsum_proposed(idx) = gbestVal;
end
% save(strcat('.\res_P\proposed.mat'), 'Pmax', 'totalsum_proposed');
save(strcat('.\res_P\proposed_SDR.mat'), 'Pmax', 'totalsum_proposed');
save(strcat('.\res_P\proposed_SDR_comput.mat'), 'Pmax', 'totalsum_proposed');