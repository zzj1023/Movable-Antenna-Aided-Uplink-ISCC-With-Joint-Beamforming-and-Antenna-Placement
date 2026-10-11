clc; 
clear; 
close all;

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
% redius = 50 + 50*rand(1,K);
% 
%  angle = 2*pi*rand(1,K);
% r_terminal = BSPos + repmat(redius,2,1).*[cos(angle);sin(angle)];
% random_param.L = 5;
% 
% random_param.theta = pi*rand(K,random_param.L) ;
% random_param.phi   = 2*pi*rand(K,random_param.L);
% random_param.phi_tx = 2*pi*rand(K,random_param.L);
% 
% random_param.Sigma = (randn(random_param.L,K)+1j*randn(random_param.L,K))/sqrt(2);
% 
% 
% % sensing target angle（一次随机，以后固定）
% random_param.theta_target = pi/6 ;
% 
% % interference random channels（一次固定）
% for k = 1:K
%     for j = 1:K
%         if k~=j
%             random_param.H_int{k,j} = (randn(N) + 1j*randn(N))/sqrt(2);
%         else
%             random_param.H_int{k,j} = zeros(N);
%         end
%     end
% end

% ================= 位置初始化 =================
% area_len = 10 * lambda;  % 天线阵列覆盖范围

% -------- 均匀排列方案（正方形阵列） --------
% 假设 M 是完全平方数，sqrt(M) x sqrt(M) 阵列
n_side = ceil(sqrt(M));
[x_grid, y_grid] = meshgrid(linspace(-area_len/2, area_len/2, n_side), linspace(-area_len/2, area_len/2, n_side));
ant_pos = [x_grid(:), y_grid(:)];
ant_pos = ant_pos(1:M,:).';   % 保证维度是 2 x M


% 设备天线线性排列
% d_tx = (0:N-1).' * (lambda/2);
% 
% save('.\channel\common_scenario_K6.mat', 'angle','r_terminal', 'd_tx', 'random_param');

%% ================= 迭代次数设置 =================
maxIter = 50;
totalsum_FPA = zeros(maxIter,1);  % Uniform Rectangular Array

%% ================= 主循环: 均匀基站天线方案 =================
for I =  1:maxIter
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

    % 调用 FP+WMMSE+SCA 核心函数

    [sum_rate, rates, SINRs, Vcell] = run_SDR_sca_sum_rate_ISAC( ...
    ant_pos, user_pos, tx_pos_rel, K, N, M, lambda, PL0, alpha_pathloss, ...
    noise_pow, Pmax,random_param, BSPos, L);

    totalsum_FPA(I) = sum_rate;

    fprintf('FPA 固定天线总速率 = %.4f bps/Hz \n', sum_rate);
end


%% ================= 保存结果 =================

% save(strcat('.\sum_conve\FPA.mat'), 'totalsum_FPA');
save(strcat('.\sum_conve\FPA_comput.mat'), 'totalsum_FPA');