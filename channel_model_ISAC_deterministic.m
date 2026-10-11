function [H, A, H_int] = channel_model_ISAC_deterministic( ...
    user_pos, ant_pos, tx_pos_rel, lambda, PL0, alpha_pathloss, ...
    random_param, K, BSPos,L)

% L = random_param.L;
M = size(ant_pos,2);
N = length(tx_pos_rel);

H = zeros(M,N,K);
H_int = cell(K,K);
d0 = 1;                      % reference 1 m
g0_lin = 10^(PL0/10);
%% ===== 1. Offloading link (H) =====
for k = 1:K
    Fk = zeros(L, M);
    for m = 1:M
        xm = ant_pos(1,m); 
        ym = ant_pos(2,m);

        for l = 1:L
            % deterministic phase (从 random_param 取)
            theta_kml = random_param.theta(k,l);
            phi_kml   = random_param.phi(k,l);

            rho_r = xm*sin(theta_kml)*cos(phi_kml) + ym*cos(theta_kml);
            Fk(l,m) = exp(1j*2*pi/lambda * rho_r);
        end
    end
    Gk = zeros(L, N);

    for n = 1:N
        xkn = tx_pos_rel(n);

        for l = 1:L
            phi_knl = random_param.phi_tx(k,l);
            rho_t = xkn * cos(phi_knl);
            Gk(l,n) = exp(1j*2*pi/lambda * rho_t);
        end
    end

    % Pathloss deterministic
    d = norm(user_pos(k,:) - BSPos.');
    % PL_dB = PL0 - 10*alpha_pathloss*log10(d/d0);
    % PL = 10^(PL_dB/10);
    % PL = PL0 * d^(-alpha_pathloss);
    % Sigma_k = diag( sqrt(PL/L) * random_param.Sigma(k) );

    PL_linear = g0_lin * (d/d0)^(-alpha_pathloss);  % g0 * d^{-alpha}
    sigmas_k = random_param.Sigma(:, k); 
    Sigma_k = diag( sqrt(PL_linear/L) * sigmas_k );
    % Sigma_k = diag( sqrt(PL_linear/L) * random_param.Sigma(k) );


    H(:,:,k) = (Gk' * Sigma_k * Fk)';
end

%% ===== 2. Sensing matrix (A) =====
theta_target = random_param.theta_target;
a_theta = exp(1j*2*pi*0.5*(0:N-1)'*sin(theta_target));
A = a_theta * a_theta';

%% ===== 3. Interference channels =====
for k = 1:K
    for j = 1:K
        if k==j
            H_int{k,j} = zeros(N,N);
        else
            dk_j = norm(user_pos(k,:) - user_pos(j,:));
            rho0 = 10^(PL0/10);
            PL_linear = rho0 / (dk_j^2);
            if isfield(random_param, 'H_int') && ...
                    size(random_param.H_int, 1) >= k && ...
                    size(random_param.H_int, 2) >= j
                H_tmp = random_param.H_int{k, j};
            else
                stream = RandStream('mt19937ar', 'Seed', k * 100 + j);
                H_tmp = (randn(stream, N, N) + ...
                    1j * randn(stream, N, N)) / sqrt(2);
            end

            H_int{k,j} = sqrt(PL_linear) * H_tmp;
        end
    end
end

end
