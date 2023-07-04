function [W_x, T, E_Omg] = h_conceFT(x, be, gam, frange, voice, fs, J, N)
    %%
    oct = log2(frange(end)/frange(1));
    F = frange(1).*2.^((0:(oct*voice))'./voice);
    f_c = (be/gam)^(1/gam)/(2*pi);
    s = f_c./(F./fs);
    
    %%
    r = rand(N, J);
    r = r./sum(r, 2);

    %%
    try
        h_fun = @(N) abs(trapz(s(1).*linspace(0, 1, N), morsewavelet(gam, be, 0, s(1).*linspace(0, 1, N)).^2) - 1) - 1e-4;
        L = round(fzero(h_fun, [numel(x)*0.5, fs]));
    catch
        L = numel(x);
    end

    ZEROPAD = zeros(1, ceil((L - numel(x))/2));
    z = [ZEROPAD, x, ZEROPAD];
    f = linspace(0, 1, numel(z));

    %% vairables
    W = zeros(numel(s), numel(z));
    dW = zeros(numel(s), numel(z));
    T = zeros(numel(s), numel(x));
    Z = fft(z);
    H_ini = zeros(1, numel(f));
    k = 0:(J - 1);

    %% GPU
    %{
    W = gpuArray(cast(W, 'single'));
    dW = gpuArray(cast(dW, 'single'));
    T = gpuArray(cast(T, 'single'));
    Z = gpuArray(cast(Z, 'single'));

    f = gpuArray(cast(f, 'single'));
    s = gpuArray(cast(s, 'single'));
    r = gpuArray(cast(r, 'single'));
    fs = gpuArray(cast(fs, 'single'));
    H_ini = gpuArray(cast(H_ini, 'single'));
    k = gpuArray(cast(k, 'single'));
    %}

    %%
    E_Omg = 0;
    fw = waitbar(0, sprintf('%s: Processing...', datetime));
    
    for n=1:N
        waitbar(n/N, fw, sprintf('%s: Processing...', datetime));

        %% orthonormal wavelets
        for i=1:numel(s)
            H = H_ini;
            dH = H_ini;
            f_s = s(i).*f;
            for j=1:J
                [H_j, dH_j] = morsewavelet(gam, be, k(j), f_s);
                H = H + r(n, j).*H_j;
                dH = dH + r(n, j).*dH_j;
            end
    
            W(i, :) = ifft(Z.*H);
            dW(i, :) = ifft(Z.*dH)./s(i);
        end
        
        %% IF estimation
        Omg = dW./W;
        Omg = Omg(:, numel(ZEROPAD) + 1:end - numel(ZEROPAD));
        Omg = fs.*real(Omg);

        E_Omg = E_Omg + Omg;
    
        %% thresholding
        W_x = W(:, numel(ZEROPAD) + 1:end - numel(ZEROPAD));
        W_x = W_x.*sqrt(s);
        P = abs(W_x);
        dlt = median(P(:))/0.6245*2;

        %% synchrosqueezing
        for i=2:(size(T, 1) - 1)
            for j=1:size(T, 2)
                if P(i, j) > dlt
                    idx = voice*log2(Omg(i, j)/frange(1)) + 1;
                    idx = round(idx);
                    if Omg(i, j) > 0 && idx > 1 && idx < numel(s)
                        c_idx = (s(idx + 1) - s(idx - 1))/s(idx)^1.5;
                        c_i = (s(i + 1) - s(i - 1))/s(i)^1.5;
                        T(idx, j) = T(idx, j) + W_x(i, j) * c_i/c_idx;
                    end
                end
            end
        end
    end
    
    T = T./N;
    E_Omg = E_Omg./N;

    close(fw);

    %% plain wavelet transform
    k = 0;
    for i=1:numel(s)
        H = morsewavelet(gam, be, k, s(i).*f);
        W(i, :) = ifft(Z.*H);
    end

    W_x = W(:, numel(ZEROPAD) + 1:end - numel(ZEROPAD));
    W_x = W_x.*sqrt(s);
end