clc;
clear;
close all;

% GA benchmark: only the outer MA-position optimizer differs from main_ISCC.
rng(1);
load('.\channel\common_scenario_K4.mat');

%% ================= System parameters =================
K = 4;
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

user_pos = r_terminal.';
tx_pos_rel = d_tx;

%% ================= GA parameters =================
ga.numGenerations = 50;       % Same outer iteration count as PSO
ga.populationSize = 70;       % Same population size as PSO
ga.dim = 2 * M;
ga.crossoverRate = 0.9;
ga.mutationRate = 1 / ga.dim;
ga.eliteCount = 2;
ga.tournamentSize = 3;
ga.blendAlpha = 0.5;
ga.sigmaMax = 0.10 * area_len;
ga.sigmaMin = 0.01 * area_len;
ga.lowerBound = -area_len / 2;
ga.upperBound = area_len / 2;

%% ================= Feasible population =================
population = zeros(ga.populationSize, ga.dim);

% Use the same type of initial point as the PSO implementation.
initial_pos = generate_MA_positions(M, lambda, area_len);
population(1, :) = initial_pos(:).';

for p = 2:ga.populationSize
    ant_pos = generate_random_MA_positions(M, area_len, d_min);
    population(p, :) = ant_pos(:).';
end

fitness = -inf(ga.populationSize, 1);
best_rate_GA = -inf;
best_pos_GA = initial_pos;
totalsum_GA = zeros(ga.numGenerations, 1);

%% ================= GA optimization =================
disp('================= Start GA benchmark =================');
disp('Outer layer: GA optimizes MA positions');
disp('Inner layer: the unchanged SDR-SCA module optimizes beamforming');

for generation = 1:ga.numGenerations
    % Every chromosome is evaluated by the same inner solver used by PSO.
    for p = 1:ga.populationSize
        ant_pos = reshape(population(p, :), [2, M]);
        fitness(p) = run_SDR_sca_sum_rate_ISAC( ...
            ant_pos, user_pos, tx_pos_rel, K, N, M, lambda, PL0, ...
            alpha_pathloss, noise_pow, Pmax, random_param, BSPos, L);
    end

    [generation_best, best_idx] = max(fitness);
    if generation_best > best_rate_GA
        best_rate_GA = generation_best;
        best_pos_GA = reshape(population(best_idx, :), [2, M]);
    end
    totalsum_GA(generation) = best_rate_GA;

    fprintf('GA generation %d/%d: best sum-rate = %.4f bps/Hz\n', ...
        generation, ga.numGenerations, best_rate_GA);

    if generation == ga.numGenerations
        break;
    end

    % Elitism keeps the best chromosomes found in the current generation.
    [~, order] = sort(fitness, 'descend');
    next_population = zeros(size(population));
    next_population(1:ga.eliteCount, :) = ...
        population(order(1:ga.eliteCount), :);

    sigma = ga.sigmaMax - (ga.sigmaMax - ga.sigmaMin) * ...
        (generation - 1) / max(ga.numGenerations - 1, 1);

    next_idx = ga.eliteCount + 1;
    while next_idx <= ga.populationSize
        parent1 = tournament_select(population, fitness, ga.tournamentSize);
        parent2 = tournament_select(population, fitness, ga.tournamentSize);

        if rand <= ga.crossoverRate
            [child1, child2] = blend_crossover( ...
                parent1, parent2, ga.blendAlpha);
        else
            child1 = parent1;
            child2 = parent2;
        end

        child1 = gaussian_mutation(child1, ga.mutationRate, sigma);
        child2 = gaussian_mutation(child2, ga.mutationRate, sigma);

        child1 = repair_ma_chromosome( ...
            child1, M, area_len, d_min, ga.lowerBound, ga.upperBound);
        child2 = repair_ma_chromosome( ...
            child2, M, area_len, d_min, ga.lowerBound, ga.upperBound);

        next_population(next_idx, :) = child1;
        if next_idx + 1 <= ga.populationSize
            next_population(next_idx + 1, :) = child2;
        end
        next_idx = next_idx + 2;
    end

    population = next_population;
end

if ~exist('.\sum_conve', 'dir')
    mkdir('.\sum_conve');
end
save('.\sum_conve\GA_1.mat', ...
    'totalsum_GA', 'best_rate_GA', 'best_pos_GA', 'ga');

%% ================= Plot =================
figure;
plot(1:ga.numGenerations, totalsum_GA, '-sg', ...
    'LineWidth', 1.5, 'MarkerFaceColor', 'w');
xlabel('GA generation');
ylabel('Sum-rate (bps/Hz)');
title('GA-based MA Position Optimization');
grid on;

%% ================= Local functions =================
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
mutation_mask = rand(size(chromosome)) <= mutation_rate;
chromosome(mutation_mask) = chromosome(mutation_mask) + ...
    sigma * randn(1, nnz(mutation_mask));
end

function chromosome = repair_ma_chromosome( ...
    chromosome, M, area_len, d_min, lower_bound, upper_bound)
chromosome = min(max(chromosome, lower_bound), upper_bound);
ant_pos = reshape(chromosome, [2, M]);

% Keep each previously repaired antenna fixed and relocate only antennas
% that violate the minimum-distance constraint.
max_attempts = 500;
for m = 2:M
    distances = vecnorm(ant_pos(:, 1:m-1) - ant_pos(:, m), 2, 1);
    attempt = 0;
    while any(distances < d_min) && attempt < max_attempts
        ant_pos(:, m) = area_len * (rand(2, 1) - 0.5);
        distances = vecnorm(ant_pos(:, 1:m-1) - ant_pos(:, m), 2, 1);
        attempt = attempt + 1;
    end

    if any(distances < d_min)
        ant_pos = generate_random_MA_positions(M, area_len, d_min);
        break;
    end
end

chromosome = ant_pos(:).';
end
