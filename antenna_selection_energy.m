function [sel_idx, sel_pos, scores] = antenna_selection_energy( ...
    ant_pos_candidate, m_active, user_pos, tx_pos_rel, ...
    lambda, PL0, alpha_pathloss, random_param, K, N, BSPos, L)

% ----------------------------------------------------
% Energy-based Antenna Pre-selection
% ----------------------------------------------------

M_all = size(ant_pos_candidate, 2);
scores = zeros(M_all, 1);

for m = 1:M_all
    r_MA = ant_pos_candidate(:, m); % 单个天线位置 2x1
    
    % 调用信道模型
    % 注意：H 的维度是 [1, N, K] (1根天线, N个发射, K个用户)
    [H, ~, ~] = channel_model_ISAC_deterministic( ...
        user_pos, r_MA, tx_pos_rel, lambda, PL0, alpha_pathloss, ...
        random_param, K, BSPos, L);

    energy_sum = 0;
    % 对所有 K 用户累积能量
    for k = 1:K
        % 提取第 k 个用户的信道 (1 x N)
        H_k = H(1, :, k); 
        % 计算该位置对用户 k 的贡献能量
        energy_sum = energy_sum + norm(H_k, 'fro')^2;
    end

    scores(m) = energy_sum;
end

% 返回得分最高的前 m_active 个位置
[~, idx_sort] = sort(scores, 'descend');
sel_idx = idx_sort(1:m_active);
sel_pos = ant_pos_candidate(:, sel_idx);

end