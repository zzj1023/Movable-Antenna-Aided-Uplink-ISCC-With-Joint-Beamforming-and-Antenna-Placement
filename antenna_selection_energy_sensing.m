function [sel_idx, sel_pos, scores] = antenna_selection_energy_sensing( ...
    ant_pos_candidate, m_active, user_pos, tx_pos_rel, ...
    lambda, PL0, alpha_pathloss, random_param, K, N,BSPos)
% ----------------------------------------------------
% Energy-based Antenna Pre-selection
% 根据总信道能量对天线进行评分
% ----------------------------------------------------

M_all = size(ant_pos_candidate, 2);
scores = zeros(M_all,1);

for m = 1:M_all

    r_MA = ant_pos_candidate(:,m);   % 单个天线位置 2×1

    energy_sum = 0;
    eta = 0.1; 

    % 对所有 K 用户累积能量
    for k = 1:K
        % H = calc_channel_singleMA( ...
        %         r_MA, user_pos(k,:).', tx_pos_rel, ...
        %         lambda, PL0, alpha_pathloss, random_param, k, N,BSPos);
        [H, A, H_int] = channel_model_ISAC_deterministic( ...
    user_pos, r_MA, tx_pos_rel, lambda, PL0, alpha_pathloss, ...
    random_param, K, BSPos);
        energy_comm  = norm(H,'fro')^2;
        energy_sense = real(trace(A' * A));
        energy_sum = energy_sum + energy_comm + eta * energy_sense;

        % energy_sum = energy_sum + norm(H,'fro')^2;
    end

    scores(m) = energy_sum;
end

% 返回前 m_active 个（但 greedy 第二步会再处理）
[~, idx_sort] = sort(scores,'descend');
sel_idx = idx_sort(1:m_active);
sel_pos = ant_pos_candidate(:, sel_idx);

end
