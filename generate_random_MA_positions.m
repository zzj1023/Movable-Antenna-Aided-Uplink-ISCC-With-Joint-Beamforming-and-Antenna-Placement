function ant_pos = generate_random_MA_positions(M, area_len, d_min)
% Generate random movable antenna positions in [-area_len/2, area_len/2]^2.
% A simple rejection sampler is used to keep pairwise distance >= d_min.

maxTry = 5000;
ant_pos = zeros(2, M);

for m = 1:M
    valid = false;
    tryCnt = 0;

    while ~valid && tryCnt < maxTry
        tryCnt = tryCnt + 1;
        newPos = area_len * (rand(2, 1) - 0.5);

        if m == 1
            valid = true;
        else
            dist_to_existing = vecnorm(ant_pos(:, 1:m-1) - newPos, 2, 1);
            valid = all(dist_to_existing >= d_min);
        end
    end

    if ~valid
        warning('Random MA position generation relaxed at antenna %d.', m);
        newPos = area_len * (rand(2, 1) - 0.5);
    end

    ant_pos(:, m) = newPos;
end
end
