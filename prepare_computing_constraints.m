function [computing, gamma_latency] = prepare_computing_constraints( ...
    computing, K)
%PREPARE_COMPUTING_CONSTRAINTS Validate CPU allocations and derive SINR limits.

required_fields = {'B_u', 'Lk', 'Ck', 'fe', 'f_total', 'Tmax'};
for idx = 1:numel(required_fields)
    field_name = required_fields{idx};
    if ~isfield(computing, field_name)
        error('Computing configuration is missing field "%s".', field_name);
    end
end

vector_fields = {'Lk', 'Ck', 'fe', 'Tmax'};
for idx = 1:numel(vector_fields)
    field_name = vector_fields{idx};
    value = computing.(field_name);
    if isscalar(value)
        value = repmat(value, K, 1);
    end
    value = value(:);
    if numel(value) ~= K || any(~isfinite(value))
        error('computing.%s must contain one finite value per device.', ...
            field_name);
    end
    computing.(field_name) = value;
end

if ~isscalar(computing.B_u) || ~isfinite(computing.B_u) || computing.B_u <= 0
    error('computing.B_u must be a positive finite scalar.');
end
if ~isscalar(computing.f_total) || ...
        ~isfinite(computing.f_total) || computing.f_total <= 0
    error('computing.f_total must be a positive finite scalar.');
end
if any(computing.Lk <= 0) || any(computing.Ck < 0) || ...
        any(computing.fe <= 0) || any(computing.Tmax <= 0)
    error('Lk, fe, and Tmax must be positive, and Ck must be nonnegative.');
end
if sum(computing.fe) > computing.f_total * (1 + 1e-12)
    error(['Infeasible CPU allocation: sum(fe)=%.4g exceeds ', ...
        'f_total=%.4g cycles/s.'], sum(computing.fe), computing.f_total);
end

execution_time = computing.Ck ./ computing.fe;
transmission_budget = computing.Tmax - execution_time;
if any(transmission_budget <= 0)
    bad_k = find(transmission_budget <= 0, 1);
    error(['Infeasible latency setting for device %d: execution time ', ...
        '%.4g s is not below Tmax %.4g s.'], bad_k, ...
        execution_time(bad_k), computing.Tmax(bad_k));
end

min_rate = computing.Lk ./ (computing.B_u * transmission_budget);
gamma_latency = 2 .^ min_rate - 1;
if any(~isfinite(gamma_latency))
    error('The latency setting produces a nonfinite minimum SINR.');
end

end
