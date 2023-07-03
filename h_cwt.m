function W_x = h_cwt(x, be, gam, frange, voice, fs)
    %%
    oct = log2(frange(end)/frange(1));
    F = frange(1).*2.^((0:(oct*voice))'./voice);
    f_c = (be/gam)^(1/gam)/(2*pi);
    s = f_c./(F./fs);
    
    %%
    %h_fun = @(N) abs(trapz(s(1).*linspace(0, 1, N), morsewavelet(gam, be, 0, s(1).*linspace(0, 1, N)).^2) - 1) - 1e-4;
    %L = round(fzero(h_fun, [numel(x), fs]));
    L = numel(x);

    if L - numel(x) > 0
        ZEROPAD = zeros(1, ceil((L - numel(x))/2));
    else
        ZEROPAD = [];
    end
    z = [ZEROPAD, x, ZEROPAD];
    f = linspace(0, 1, numel(z));

    %%
    Z = fft(z);
    W = zeros(numel(s), numel(z));

    %% plain wavelet transform
    k = 0;
    for i=1:numel(s)
        H = morsewavelet(gam, be, k, s(i).*f);
        W(i, :) = ifft(Z.*H);
    end

    W_x = W(:, numel(ZEROPAD) + 1:end - numel(ZEROPAD));
    W_x = W_x.*sqrt(s);
end