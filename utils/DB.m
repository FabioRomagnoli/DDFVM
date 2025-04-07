function [res] = DB(x)
    [bp, bn] = bimu_bernoulli(x);

    res = zeros(size(x));

    % Handle x ≠ 0 normally
    nonzero_idx = (x ~= 0);
    res(nonzero_idx) = (bp(nonzero_idx) ./ x(nonzero_idx)) .* (1 - bn(nonzero_idx));

    % Handle x = 0 specially
    zero_idx = (x == 0);
    res(zero_idx) = -0.5;
end
