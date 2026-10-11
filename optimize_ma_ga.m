function [best_rate, best_pos, best_history] = optimize_ma_ga( ...
    user_pos, tx_pos_rel, K, N, M, lambda, PL0, alpha_pathloss, ...
    noise_pow, Pmax, random_param, BSPos, L, area_len, d_min, ...
    ga, computing, initial_pos)
%OPTIMIZE_MA_GA Optimize MA positions while keeping SDR-SCA unchanged.

dim = 2 * M;
lower_bound = -area_len / 2;
upper_bound = area_len / 2;

population = zeros(ga.populationSize, dim);
if isempty(initial_pos)
    initial_pos = generate_MA_positions(M, lambda, area_len);
end
population(1, :) = initial_pos(:).';
for p = 2:ga.populationSize
    candidate = generate_random_MA_positions(M, area_len, d_min);
    population(p, :) = candidate(:).';
end

fitness = -inf(ga.populationSize, 1);
best_rate = -inf;
best_pos = initial_pos;
best_history = zeros(ga.numGenerations, 1);

for generation = 1:ga.numGenerations
    for p = 1:ga.populationSize
        ant_pos = reshape(population(p, :), [2, M]);
        fitness(p) = run_SDR_sca_sum_rate_ISAC( ...
            ant_pos, user_pos, tx_pos_rel, K, N, M, lambda, PL0, ...
            alpha_pathloss, noise_pow, Pmax, random_param, BSPos, L, ...
            computing);
    end

    [generation_best, best_idx] = max(fitness);
    if generation_best > best_rate
        best_rate = generation_best;
        best_pos = reshape(population(best_idx, :), [2, M]);
    end
    best_history(generation) = best_rate;

    fprintf('GA generation %d/%d | Best sum-rate = %.4f bps/Hz\n', ...
        generation, ga.numGenerations, best_rate);

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
        child1 = repair_position( ...
            child1, M, area_len, d_min, lower_bound, upper_bound);
        child2 = repair_position( ...
            child2, M, area_len, d_min, lower_bound, upper_bound);

        next_population(next_idx, :) = child1;
        if next_idx + 1 <= ga.populationSize
            next_population(next_idx + 1, :) = child2;
        end
        next_idx = next_idx + 2;
    end
    population = next_population;
end

end

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
    rand(size(parent2)) .* ((1 + 2 * alpha) * span);
end

function child = gaussian_mutation(child, mutation_rate, sigma)
mutation_mask = rand(size(child)) < mutation_rate;
child(mutation_mask) = child(mutation_mask) + ...
    sigma * randn(1, nnz(mutation_mask));
end

function chromosome = repair_position( ...
    chromosome, M, area_len, d_min, lower_bound, upper_bound)
chromosome = min(max(chromosome, lower_bound), upper_bound);
ant_pos = reshape(chromosome, [2, M]);

for m = 2:M
    attempts = 0;
    while any(vecnorm(ant_pos(:, 1:m-1) - ant_pos(:, m), 2, 1) < d_min)
        ant_pos(:, m) = lower_bound + area_len * rand(2, 1);
        attempts = attempts + 1;
        if attempts >= 200
            ant_pos = generate_random_MA_positions(M, area_len, d_min);
            chromosome = ant_pos(:).';
            return;
        end
    end
end
chromosome = ant_pos(:).';
end
