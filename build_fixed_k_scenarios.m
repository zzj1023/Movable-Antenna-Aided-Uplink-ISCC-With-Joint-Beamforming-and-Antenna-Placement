clc;
clear;

source_file = fullfile('.\scenario_saved', ...
    'fixed_users_and_channels.mat');
scenario = load(source_file);

Ks = [2, 4, 6, 8, 10];
N = 4;
M = 4;
lambda = 0.1;

if scenario.Kmax < max(Ks)
    error('The source scenario supports at most K=%d.', scenario.Kmax);
end

% All K values use the same transmit array and fixed receive array.
d_tx = (0:N-1).' * (lambda / 2);
d = lambda / 2;
n_side = ceil(sqrt(M));
[xg, yg] = meshgrid(0:n_side-1, 0:n_side-1);
xg = xg(:) * d;
yg = yg(:) * d;
ant_pos = [xg(1:M), yg(1:M)].';
ant_pos = ant_pos - mean(ant_pos, 2);

for idx = 1:numel(Ks)
    K = Ks(idx);

    % Nested scenarios: the K-user case always uses the first K users
    % and the corresponding prefix of the same Kmax-user channel bank.
    user_pos = scenario.user_pos_all(1:K, :);
    r_terminal = scenario.r_terminal_all(:, 1:K);

    random_param.L = scenario.random_big.L;
    random_param.theta = scenario.random_big.theta(1:K, :);
    random_param.phi = scenario.random_big.phi(1:K, :);
    random_param.phi_tx = scenario.random_big.phi_tx(1:K, :);
    random_param.Sigma = scenario.random_big.Sigma(:, 1:K);
    random_param.theta_target = scenario.random_big.theta_target;
    random_param.H_int = scenario.random_big.H_int(1:K, 1:K);

    computing = default_computing_config(K);
    Dk = computing.Lk;
    Ck = computing.Ck;
    fe_k = computing.fe;

    save_dir = fullfile('.\scenario_saved', ['K', num2str(K)]);
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    save(fullfile(save_dir, 'fixed_env.mat'), ...
        'K', 'user_pos', 'r_terminal', 'Dk', 'Ck', 'fe_k', 'computing', ...
        'd_tx', 'random_param', 'ant_pos');

    fprintf('Saved nested fixed environment for K=%d.\n', K);
end
