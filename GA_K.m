clc;
clear;
close all;
rng(1);

%% System parameters
Ks = [2, 4, 6, 8, 10];
N = 4;
M = 4;
L = 5;
lambda = 0.1;
area_len = 8 * lambda;
d_min = lambda / 2;
noise_pow = 1e-9;
Pmax = 1;
PL0 = -40;
alpha_pathloss = 2.8;
BSPos = [0; 0];

%% GA parameters: same settings as main_GA
ga.numGenerations = 50;
ga.populationSize = 70;
ga.dim = 2 * M;
ga.crossoverRate = 0.9;
ga.mutationRate = 1 / ga.dim;
ga.eliteCount = 2;
ga.tournamentSize = 3;
ga.blendAlpha = 0.5;
ga.sigmaMax = 0.10 * area_len;
ga.sigmaMin = 0.01 * area_len;

GA_sumRate = zeros(numel(Ks), 1);

for idx = 1:numel(Ks)
    K = Ks(idx);
    env_file = fullfile('.\scenario_saved', ...
        ['K', num2str(K)], 'fixed_env.mat');
    env = load(env_file);
    user_pos = env.user_pos;
    tx_pos_rel = env.d_tx;
    current_param = env.random_param;
    computing = default_computing_config(K);
    computing.Lk = env.Dk;
    computing.Ck = env.Ck;
    computing.fe = env.fe_k;

    % Same MA initialization used by main_GA.
    population = zeros(ga.populationSize, ga.dim);
    initial_pos = generate_MA_positions(M, lambda, area_len);
    population(1, :) = initial_pos(:).';
    for p = 2:ga.populationSize
        candidate = generate_random_MA_positions(M, area_len, d_min);
        population(p, :) = candidate(:).';
    end

    fitness = -inf(ga.populationSize, 1);
    best_rate = -inf;

    for generation = 1:ga.numGenerations
        for p = 1:ga.populationSize
            ant_pos = reshape(population(p, :), [2, M]);
            fitness(p) = run_SDR_sca_sum_rate_ISAC( ...
                ant_pos, user_pos, tx_pos_rel, K, N, M, lambda, ...
                PL0, alpha_pathloss, noise_pow, Pmax, current_param, ...
                BSPos, L, computing);
        end

        generation_best = max(fitness);
        best_rate = max(best_rate, generation_best);
        fprintf('K=%d | GA generation=%d/%d | Best sum-rate=%.4f bps/Hz\n', ...
            K, generation, ga.numGenerations, best_rate);

        if generation == ga.numGenerations
            break;
        end

        [~, order] = sort(fitness, 'descend');
        next_population = zeros(size(population));
        next_population(1:ga.eliteCount, :) = ...
            population(order(1:ga.eliteCount), :);

        sigma = ga.sigmaMax - (ga.sigmaMax - ga.sigmaMin) * ...
            (generation - 1) / max(ga.numGenerations - 1, 1);
        next_idx = ga.eliteCount + 1;

        while next_idx <= ga.populationSize
            parent1 = tournament_select( ...
                population, fitness, ga.tournamentSize);
            parent2 = tournament_select( ...
                population, fitness, ga.tournamentSize);

            if rand <= ga.crossoverRate
                [child1, child2] = blend_crossover( ...
                    parent1, parent2, ga.blendAlpha);
            else
                child1 = parent1;
                child2 = parent2;
            end

            child1 = gaussian_mutation( ...
                child1, ga.mutationRate, sigma);
            child2 = gaussian_mutation( ...
                child2, ga.mutationRate, sigma);
            child1 = repair_ma_position(child1, M, area_len, d_min);
            child2 = repair_ma_position(child2, M, area_len, d_min);

            next_population(next_idx, :) = child1;
            if next_idx + 1 <= ga.populationSize
                next_population(next_idx + 1, :) = child2;
            end
            next_idx = next_idx + 2;
        end

        population = next_population;
    end

    GA_sumRate(idx) = best_rate;
end

if ~exist('.\res_K', 'dir')
    mkdir('.\res_K');
end
save('.\res_K\GA_SDR.mat', 'Ks', 'GA_sumRate', 'ga');

function parent = tournament_select(population, fitness, tournament_size)
candidate_idx = randi(size(population, 1), tournament_size, 1);
[~, local_best] = max(fitness(candidate_idx));
parent = population(candidate_idx(local_best), :);
end

function [child1, child2] = blend_crossover(parent1, parent2, alpha)
lower = min(parent1, parent2);
upper = max(parent1, parent2);
span = upper - lower;
child1 = lower - alpha * span + ...
    rand(size(parent1)) .* ((1 + 2 * alpha) * span);
child2 = lower - alpha * span + ...
    rand(size(parent1)) .* ((1 + 2 * alpha) * span);
end

function chromosome = gaussian_mutation(chromosome, mutation_rate, sigma)
mask = rand(size(chromosome)) <= mutation_rate;
chromosome(mask) = chromosome(mask) + sigma * randn(1, nnz(mask));
end

function chromosome = repair_ma_position(chromosome, M, area_len, d_min)
ant_pos = reshape(chromosome, [2, M]);
ant_pos = min(max(ant_pos, -area_len / 2), area_len / 2);

for m = 2:M
    distances = vecnorm(ant_pos(:, 1:m-1) - ant_pos(:, m), 2, 1);
    attempts = 0;
    while any(distances < d_min) && attempts < 500
        ant_pos(:, m) = area_len * (rand(2, 1) - 0.5);
        distances = vecnorm(ant_pos(:, 1:m-1) - ant_pos(:, m), 2, 1);
        attempts = attempts + 1;
    end
end

chromosome = ant_pos(:).';
end
