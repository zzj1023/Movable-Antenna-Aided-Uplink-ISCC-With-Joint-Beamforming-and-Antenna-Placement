function [r_MA] = generate_MA_positions(M, lambda, area_len)

minDist = lambda/2;     % 任意两天线最小距离要求
r_MA = zeros(2, M); 

for m = 1:M
    valid = false;
    while ~valid
        % 随机生成新点
        newPos = area_len * (rand(2,1) - 0.5);

        if m == 1
            % 第一个点直接放
            r_MA(:,m) = newPos;
            valid = true;
        else
            % 检查与所有已放置点的距离
            dist_to_existing = vecnorm(r_MA(:,1:m-1) - newPos, 2, 1);
            if all(dist_to_existing >= minDist)
                r_MA(:,m) = newPos;
                valid = true;
            end
        end
    end
end
end
