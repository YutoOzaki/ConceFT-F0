function [W_x, T, E_Omg] = h_conceFT2(x, be, gam, frange, voice, fs, J, N)
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
    Omg = zeros(numel(s), numel(z));
    T = zeros(numel(s), numel(x));
    X = fft(z);
    H_ini = zeros(1, numel(f));
    k = 0:(J - 1);
    t = 0:(numel(x) - 1);

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
            xiH = H_ini;
            xisqH = H_ini;
            dH = H_ini;
            xidH = H_ini;

            for j=1:J
                [H_j, xiH_j, xisqH_j, dH_j, xidH_j] = morsewavelet(gam, be, k(j), s(i).*f);
                
                H = H + r(n, j).*H_j;
                xiH = xiH + r(n, j).*xiH_j;
                xisqH = xisqH + r(n, j).*xisqH_j;
                dH = dH + r(n, j).*dH_j;
                xidH = xidH + r(n, j).*xidH_j;
            end
    
            W_H = ifft(X.*H).*sqrt(s(i));
            W_xiH = ifft(X.*xiH).*sqrt(s(i));
            W_xisqH = ifft(X.*xisqH).*sqrt(s(i));
            W_dH = ifft(X.*dH).*sqrt(s(i));
            W_xidH = ifft(X.*xidH).*sqrt(s(i));

            W(i, :) = W_H;
            
            Omg_1 = (1/s(i)).*W_xiH./W_H;
            tau = t + s(i)/(1i*2*pi).*(W_dH./W_H);

            dOmg = (1i*2*pi/s(i)) .* (W_H.*W_xisqH - W_xiH.^2)./W_H.^2;
            dtau = 1/fs + s(i).*(W_H.*W_xidH - W_dH.*W_xiH)./W_H.^2;
            q = dOmg./dtau;
            Omg_2 = Omg_1 + q.*(t - tau);
        
            Omg(i, :) = Omg_2;
        end
        
        %% IF estimation
        Omg = Omg(:, numel(ZEROPAD) + 1:end - numel(ZEROPAD));
        Omg = fs.*real(Omg);

        E_Omg = E_Omg + Omg;
    
        %% thresholding
        W_x = W(:, numel(ZEROPAD) + 1:end - numel(ZEROPAD));
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
        W(i, :) = ifft(X.*H);
    end

    W_x = W(:, numel(ZEROPAD) + 1:end - numel(ZEROPAD));
    W_x = W_x.*sqrt(s);
end