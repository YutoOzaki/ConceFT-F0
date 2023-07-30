function y = h_denoising(s, fs)
    %%
    N = 2^nextpow2(fs*0.025);
    window = sqrt(hann(N, 'periodic'));

    L = round(fs*0.01);
    J = 4096;
    q_u = 0.2;
    q_l = 0.1;
    
    %%
    E = zeros(numel(s), 1);
    for i=1:numel(s)
        idx_st = max(1, i - L);
        idx_ed = min(numel(s), i + L);
        E(i) = mean(s(idx_st:idx_ed).^2);
    end
    
    idx = find(E < quantile(E, q_u) & E > quantile(E, q_l));
    idx = idx(idx > N/2);
    idx = idx(idx < (numel(s) - N/2));
    
    A_d = zeros(N, 1);
    for j=1:J
        i = idx(randi(numel(idx)));
        x = s(i - N/2:i + N/2 - 1);
        A_d = A_d + abs(fft(window.*x));
    end
    A_d = A_d./J;
    
    %%
    noverlap = N*1/2;
    
    assert(iscola(window, noverlap), 'Check the COLA condition');
    S = stft(s, fs, 'Window', window, 'OverlapLength', noverlap, 'FrequencyRange', 'twosided');
    A_y = abs(S);
    
    y = s.*0;
    M = N - noverlap;
    K = 2*pi*(0:(N - 1))'./N.*(0:(N - 1));
    for i=1:size(A_y, 2)
        idx_st = 1 + (i - 1)*M;
        idx_ed = idx_st + N - 1;
        H = (A_y(:, i) - A_d)./A_y(:, i);
        A_s = H.*A_y(:, i);
        A_s(A_s < 0) = 0;
        A = A_s;
        theta = angle(S(:, i));
        y_i = sum(A.*cos(K + theta), 1)'./N;
        y(idx_st:idx_ed) = y(idx_st:idx_ed) + window.*y_i;
    end
    y(isnan(y)) = 0;
end