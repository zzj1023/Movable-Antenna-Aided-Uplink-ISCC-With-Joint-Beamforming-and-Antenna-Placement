function computing = default_computing_config(K)
%DEFAULT_COMPUTING_CONFIG Common edge-computing parameters for all schemes.

computing.B_u = 10e6;                    % Uplink bandwidth (Hz)
computing.Lk = 2e6 * ones(K, 1);         % Input data size (bit)
computing.vartheta = 200 * ones(K, 1);   % CPU cycles per input bit
computing.Ck = computing.vartheta .* computing.Lk;
computing.fe = 3e9 * ones(K, 1);         % Predetermined CPU allocation (cycles/s)
computing.f_total = 30e9;                % Edge-server CPU budget (cycles/s)
computing.Tmax = 5 * ones(K, 1);         % Maximum end-to-end latency (s)

end
