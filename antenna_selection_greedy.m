function [sel_idx, sel_pos, rate_hist] = antenna_selection_greedy( ...
        ant_pos_all, m_active, user_pos, tx_pos_rel, ...
        K, N, lambda, PL0, alpha_pathloss, ...
        noise_pow, Pmax, B, Dk, Ck, fe, random_param, BSPos, L)

M_all = size(ant_pos_all, 2);
avail = true(1, M_all);

computing = default_computing_config(K);
computing.B_u = B;
computing.Lk = Dk;
computing.Ck = Ck;
computing.fe = fe;

sel_idx = [];
rate_hist = zeros(m_active, 1);

for t = 1:m_active
    best_local_idx = -1;
    best_local_rate = -inf; % 初始化为极小值，我们要找最大速率

    fprintf('\n=== Greedy Step %d / %d ===\n', t, m_active);

    % 获取当前可用的候选索引
    candidates = find(avail);
    
    for cand = candidates
        % 尝试添加这根天线
        idx_try = [sel_idx, cand];
        ant_pos_try = ant_pos_all(:, idx_try); % 维度 2 x t

        % 调用评估函数 (注意传入 ant_pos_try 和当前天线数 t)
        [current_rate, ~, ~, ~] = run_SDR_sca_sum_rate_ISAC( ...
            ant_pos_try, user_pos, tx_pos_rel, K, N, t, ... 
            lambda, PL0, alpha_pathloss, noise_pow, Pmax, ...
            random_param, BSPos, L, computing);

        fprintf('  尝试索引 %d: 当前总速率 = %.4f bps/Hz\n', cand, current_rate);

        % 寻找最大速率
        if current_rate > best_local_rate
            best_local_rate = current_rate;
            best_local_idx = cand;
        end
    end

    % 更新选中的集合
    if best_local_idx ~= -1
        sel_idx = [sel_idx, best_local_idx];
        avail(best_local_idx) = false;
        rate_hist(t) = best_local_rate;
        fprintf('** 最终选中索引 %d, 此时总速率 = %.4f\n', best_local_idx, best_local_rate);
    end
end

sel_pos = ant_pos_all(:, sel_idx);

end
