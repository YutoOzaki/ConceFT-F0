function [segment, V] = h_voicingknn(y, N, minseglen)
    %%
    L = numel(y);
    E = zeros(L, 1);

    for i=1:L
        E(i) = sum(y(max(1, i - N):min(L, i + N)).^2);
    end

    %%
    V = kmeans(log10(E), 2);
    V(V == V(1)) = -1;
    V(V ~= V(1)) = 1;
    V(V == -1) = 2;

    %%
    segment = [];
    idx_ed = 0;
    idx_st = find(V(idx_ed + 1:end) == 1, 1, 'first') + idx_ed;

    while ~isempty(idx_st)
        idx_ed = find(V(idx_st:end) ~= 1, 1, 'first') + idx_st - 2;
        
        if (idx_ed - idx_st + 1) >= minseglen
            segment(end + 1, :) = [idx_st, idx_ed];
        end

        idx_st = find(V(idx_ed + 1:end) == 1, 1, 'first') + idx_ed;
    end
end