clc;
clear;
close all;

rng(1);
load('.\channel\common_scenario_K4.mat');

%% ================= System parameters =================
K = 4;
N = 4;
M = 4;
L = 5;
lambda = 0.1;
noise_pow = 1e-9;
Pmax = 1;
PL0 = -40;
alpha_pathloss = 2.8;
BSPos = [0; 0];

area_factors = [3 4 5 6 7 8];
area_list = area_factors * lambda;
FPA_sumRate = zeros(length(area_list), 1);

user_pos = r_terminal.';
if ~exist('d_tx', 'var')
    d_tx = (0:N-1).' * (lambda / 2);
end
tx_pos_rel = d_tx;

%% ================= Area sweep =================
for idx = 1:length(area_list)
    d = lambda / 2;
    n_side = ceil(sqrt(M));

    [xg, yg] = meshgrid(0:n_side-1, 0:n_side-1);
    xg = xg(:) * d;
    yg = yg(:) * d;
    ant_pos = [xg(1:M), yg(1:M)].';
    ant_pos = ant_pos - mean(ant_pos, 2);

    [FPA_sumRate(idx), ~, ~, ~] = run_SDR_sca_sum_rate_ISAC( ...
        ant_pos, user_pos, tx_pos_rel, K, N, M, lambda, PL0, alpha_pathloss, ...
        noise_pow, Pmax, random_param, BSPos, L);

    fprintf('FPA Area=%.1f lambda | Sum-Rate=%.4f bps/Hz\n', ...
        area_factors(idx), FPA_sumRate(idx));
end

if ~exist('.\res_Area', 'dir')
    mkdir('.\res_Area');
end
save(strcat('.\res_Area\FPA_SDR.mat'), 'area_factors', 'area_list', 'FPA_sumRate');
