function idx_ridge = h_viterbi(P, Dlt, c)
    %%
    R = tiedrank(-P) - 1;
    [Q, N] = size(P);
    
    %%
    C = R;
    I = R.*0;

    for n=2:N
        for q=1:Q
            idx = max(q - Dlt, 1):min(q + Dlt, Q);
            cost_zero = C(idx, n - 1);
            [cost_min, idx_min] = min(cost_zero);
            idx_cost = idx(idx_min);
            
            j = 1;
            while j*c < cost_min
                idx_l = max(q - Dlt - j, 1);
                idx_u = min(q + Dlt + j, Q);
                cost_l = C(idx_l, n - 1) + j*c;
                cost_u = C(idx_u, n - 1) + j*c;
                [cost_min, idx_min] = min([cost_min, cost_l, cost_u]);

                idx_tmp = [idx_cost, idx_l, idx_u];
                idx_cost = idx_tmp(idx_min);

                j = j + 1;
            end

            I(q, n) = idx_cost;

            if abs(idx_cost - q) <= Dlt
                penalty = 0;
            else
                penalty = c*(abs(idx_cost - q) - Dlt);
            end
            C(q, n) = C(idx_cost, n - 1) + R(q, n) + penalty;
        end
    end

    %%
    idx_ridge = zeros(N, 1);
    [~, idx_ridge(end)] = min(C(:, end));
    for n=N:-1:2
        idx_ridge(n - 1) = I(idx_ridge(n), n);
    end
end