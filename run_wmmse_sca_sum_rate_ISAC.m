function [sum_rate, rates, SINRs, Vcell, latency] = run_wmmse_sca_sum_rate_ISAC( ...
    ant_pos, user_pos, tx_pos_rel, K, N, M, lambda, PL0, alpha_pathloss, ...
    noise_pow, Pmax, random_param, BSPos, L, varargin)

% =============================================================
% ISAC beamforming module: SDR + DC-SCA
% Objective: maximize system communication sum-rate under power
% sensing, power, and edge-computing latency constraints.
%
% SDR variable:
%   W_k = v_k * v_k' , W_k >= 0, rank(W_k)=1
%
% The rank-one constraint is relaxed during CVX optimization.
% After each SDR step, the beam vector v_k is recovered from the
% dominant eigenvector of W_k.
% =============================================================

%% ---------- 0. Edge-computing parameters ----------
if isempty(varargin)
    computing = default_computing_config(K);
else
    computing = varargin{1};
end
[computing, gamma_latency] = prepare_computing_constraints(computing, K);

%% ---------- 1. Build channels ----------
[H_all, A_theta, H_int] = channel_model_ISAC_deterministic( ...
    user_pos, ant_pos, tx_pos_rel, lambda, PL0, alpha_pathloss, ...
    random_param, K, BSPos, L);

Btmp = A_theta' * A_theta;
B_sense = (Btmp + Btmp') / 2 + 1e-8 * eye(N);

%% ---------- 2. Initialize rank-one transmit covariance ----------
Vcell = cell(K, 1);
Wcell = cell(K, 1);
for k = 1:K
    heq = sum(H_all(:, :, k), 1).';
    if norm(heq) == 0 || any(~isfinite(heq))
        heq = ones(N, 1);
    end
    Vcell{k} = heq / norm(heq) * sqrt(Pmax);
    Wcell{k} = Vcell{k} * Vcell{k}';
end

Gamma_dB = 1;
Gamma = 10^(Gamma_dB / 10) * ones(K, 1);
beta_sense = 30^(-2) * ones(K, 1);
sigma2 = noise_pow * ones(K, 1);

maxOuter = 12;
tol = 1e-3;
prev_sum_rate = -inf;
stagnant_cnt = 0;
rates = zeros(K, 1);
SINRs = zeros(K, 1);
sum_rate = 0;
latency = inf(K, 1);

%% ---------- 3. SDR + DC-SCA loop ----------
for outer = 1:maxOuter
    Wprev = Wcell;

    % Step A: update linear receivers using the currently recovered beams.
    Ucell = cell(K, 1);
    for k = 1:K
        Sig = noise_pow * eye(M);
        for j = 1:K
            Sig = Sig + H_all(:, :, j) * Wprev{j} * H_all(:, :, j)';
        end

        bvec = H_all(:, :, k) * Vcell{k};
        traceSig = real(trace(Sig));
        reg_eps = max(1e-9, 1e-6 * traceSig / max(1, M));
        Sig_reg = Sig + reg_eps * eye(M);

        if rcond(Sig_reg) < 1e-12
            U_k = pinv(Sig_reg) * bvec;
        else
            U_k = Sig_reg \ bvec;
        end

        if any(~isfinite(U_k))
            U_k = (randn(M, 1) + 1j * randn(M, 1)) / sqrt(2);
        end
        Ucell{k} = U_k;
    end

    % Step B: build quadratic receive-side channel matrices.
    Q = cell(K, K);
    for rx = 1:K
        for tx = 1:K
            Q{rx, tx} = H_all(:, :, tx)' * (Ucell{rx} * Ucell{rx}') * H_all(:, :, tx);
            Q{rx, tx} = (Q{rx, tx} + Q{rx, tx}') / 2;
        end
    end

    % Step C: constants at the previous point for DC linearization.
    I0 = zeros(K, 1);
    sense_rhs = zeros(K, 1);
    for k = 1:K
        I0(k) = noise_pow * real(Ucell{k}' * Ucell{k});
        for j = 1:K
            if j ~= k
                I0(k) = I0(k) + real(trace(Q{k, j} * Wprev{j}));
            end
        end
        I0(k) = max(I0(k), 1e-12);

        sense_interference = sigma2(k);
        for j = 1:K
            if j ~= k
                Ckj = H_int{k, j}' * H_int{k, j};
                sense_interference = sense_interference + real(trace(Ckj * Wprev{j}));
            end
        end
        sense_rhs(k) = (Gamma(k) / (beta_sense(k)^2)) * sense_interference;
    end

    % Step D: solve the relaxed covariance problem.
    cvx_clear;
    cvx_begin quiet
        variable W(N, N, K) complex
        expressions total_lb total_power(K) sense_power(K)
        expressions interf(K) total_rx(K)

        total_lb = 0;
        for k = 1:K
            interf(k) = noise_pow * real(Ucell{k}' * Ucell{k});
            total_rx(k) = noise_pow * real(Ucell{k}' * Ucell{k});
            for j = 1:K
                total_rx(k) = total_rx(k) + real(trace(Q{k, j} * W(:, :, j)));
                if j ~= k
                    interf(k) = interf(k) + real(trace(Q{k, j} * W(:, :, j)));
                end
            end

            % DC lower bound:
            % log(total_rx) - log(interf)
            % >= log(total_rx) - log(I0) - (interf - I0) / I0
            total_lb = total_lb + log(total_rx(k)) / log(2) ...
                - log(I0(k)) / log(2) ...
                - (interf(k) - I0(k)) / (I0(k) * log(2));
        end

        maximize(total_lb)

        subject to
            for k = 1:K
                W(:, :, k) == hermitian_semidefinite(N);
                total_power(k) = real(trace(W(:, :, k)));
                total_power(k) <= Pmax;

                sense_power(k) = real(trace(B_sense * W(:, :, k)));
                sense_power(k) >= sense_rhs(k);

                total_rx(k) - interf(k) >= gamma_latency(k) * interf(k);
            end
    cvx_end

    if ~(strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved'))
        disp(['SDR+DC-SCA Iter ', num2str(outer), ...
              ' CVX Status: ', cvx_status, '. Keep previous beams.']);
        break;
    end

    % Step E: rank-one recovery from SDR covariance matrices.
    for k = 1:K
        Wk = full(W(:, :, k));
        Wk = (Wk + Wk') / 2;
        [Vk_eig, Dk_eig] = eig(Wk);
        [lambda_max, idx_max] = max(real(diag(Dk_eig)));
        lambda_max = max(lambda_max, 0);

        if lambda_max <= 1e-12 || any(~isfinite(Vk_eig(:, idx_max)))
            vk = Vcell{k};
        else
            vk = sqrt(lambda_max) * Vk_eig(:, idx_max);
            if norm(vk) > sqrt(Pmax)
                vk = vk / norm(vk) * sqrt(Pmax);
            end
        end

        Vcell{k} = vk;
        Wcell{k} = vk * vk';
    end

    % Step F: evaluate the feasible rank-one beams.
    [sum_rate, rates, SINRs] = evaluate_rank_one_rate( ...
        H_all, Ucell, Vcell, noise_pow, K);
    latency = computing.Lk ./ (computing.B_u * rates) + ...
        computing.Ck ./ computing.fe;

    fprintf(['SDR+DC-SCA Iter %d, Sum-Rate = %.4f bps/Hz, ', ...
        'Max Latency = %.4f s\n'], outer, sum_rate, max(latency));

    improvement = sum_rate - prev_sum_rate;
    if improvement < tol && improvement > -1e-4
        stagnant_cnt = stagnant_cnt + 1;
        if stagnant_cnt >= 2
            break;
        end
    else
        stagnant_cnt = 0;
    end
    prev_sum_rate = sum_rate;
end

end

function [sum_rate, rates, SINRs] = evaluate_rank_one_rate( ...
    H_all, Ucell, Vcell, noise_pow, K)

rates = zeros(K, 1);
SINRs = zeros(K, 1);

for k = 1:K
    num = abs(Ucell{k}' * H_all(:, :, k) * Vcell{k})^2;
    den = noise_pow * real(Ucell{k}' * Ucell{k});
    for j = 1:K
        if j ~= k
            den = den + abs(Ucell{k}' * H_all(:, :, j) * Vcell{j})^2;
        end
    end

    den = max(real(den), 1e-12);
    SINRs(k) = min(real(max(0, num / den)), 1e12);
    rates(k) = log2(1 + SINRs(k));

    if ~isfinite(rates(k)) || rates(k) < 1e-12
        rates(k) = 1e-12;
    end
end

sum_rate = sum(rates);
end
