function test07_sst2
    rng(1);

    %%
    fs = 2000;
    t = (0:(2*fs - 1))./fs;
    
    %%{
    a = rand(3, 1);
    f0 = sort((fs*0.3).*rand(3, 1));
    x = 0;
    for i=1:numel(a)
        x = x + a(i).*sin(2*pi*f0(i).*t);
    end
    %}
    
    %{
    f0 = fs*0.3*rand + 0.025*fs;
    x = chirp(t, f0, t(end), fs*0.45, 'quadratic');
    %}

    %{
    f0 = [];
    A = 1 + 3.*t.^2 + 4.*(1 - t).^7;
    ph = 240.*t - 2.*exp(-2.*t).*sin(14*pi.*t);
    x = A.*exp(2*pi*1i.*ph);
    x = real(x);
    %}

    %%
    %K = randi(5) - 1;
    K = 0;
    be = 20 + rand*10;
    gam = 2 + rand*9;
    f = linspace(0, 1, numel(x));
    X = fft(x);
    [~, ome_c] = morsefreq(be, gam);
    
    r = rand(K + 1, 1);
    r = r./sum(r);

    %%
    s_min = ome_c/0.5;
    s_max = ome_c/(1/numel(x));
    J = log2(s_max) - log2(s_min);
    dj = 1/64;
    s = s_min*2.^(0:dj:J)';

    %%
    W = zeros(numel(s), numel(x));
    Omg = zeros(numel(s), numel(x));
    t_N = linspace(0, 1, numel(x));

    for i=1:numel(s)
        H = 0;
        dH = 0;
        ddH = 0;
        kH = 0;
        dkH = 0;

        f_s = s(i).*f;

        k = 0;
        [H_k, dH_k, ddH_k] = morsewavelet(gam, be, k, f_s);
        H = H + r(k + 1).*H_k;
        dH = dH + r(k + 1).*dH_k;
        ddH = ddH + r(k + 1).*ddH_k;
        [kH_k, dkH_k, ddkH_k] = morsewavelet(gam, be, k + 1, f_s);
        kH = kH + r(k + 1).*kH_k;
        dkH = dkH + r(k + 1).*dkH_k;
        
        for k=0:K
            H_k = kH_k;
            dH_k = dkH_k;
            ddH_k = ddkH_k;

            H = H + r(k + 1).*H_k;
            dH = dH + r(k + 1).*dH_k;
            ddH = ddH + r(k + 1).*ddH_k;

            [kH_k, dkH_k, ddkH_k] = morsewavelet(gam, be, k + 1, f_s);
            kH = kH + r(k + 1).*kH_k;
            dkH = dkH + r(k + 1).*dkH_k;
        end

        W(i, :) = ifft(X.*H);
        dW = ifft(X.*dH);
        Omg_1 = dW./W(i, :)./s(i);

        kW = ifft(X.*kH);
        tau = t_N + s(i)/(2*1i*pi).*(kW./W(i, :));

        ddW = ifft(X.*ddH);
        dkW = ifft(X.*dkH);
        q = (2*1i*pi/s(i)^2) .* (ddW.*W(i, :) - dW.^2)./(W(i, :).^2 + dkW.*W(i, :) - kW.*dW);

        Omg(i, :) = Omg_1 + q.*(t_N - tau);
    end
    W = W.*sqrt(s);
    
    %% synchrosqueezing
    F = ome_c./s.*fs;
    Omg = real(Omg).*fs; % convert to frequency scale
    voice = 1/dj;
    T = zeros(numel(s), numel(x));
    
    for i=2:(size(T, 1) - 1)
        for j=1:size(T, 2)
            idx = voice*log2(F(1)/Omg(i, j)) + 1;
            idx = round(idx);
            if Omg(i, j) > 0 && idx > 1 && idx < numel(s)
                c_idx = (s(idx + 1) - s(idx - 1))/s(idx)^1.5;
                c_i = (s(i + 1) - s(i - 1))/s(i)^1.5;
                T(idx, j) = T(idx, j) + W(i, j) * c_i/c_idx;
            else
                T(i, j) = T(i, j) + W(i, j);
            end
        end
    end
    T(1, :) = T(1, :) + W(1, :);
    T(end, :) = T(end, :) + W(end, :);

    %% Reconstruction using a delta function as the synthesis wavelet
    if K == 0
        [~, C_dlt] = morseadmissibility(gam, be);
        H = morsewavelet(gam, be, K, f);
    else
        H = 0;
        for k=0:K
            H = H + r(k + 1).*morsewavelet(gam, be, k, f);
        end
        C_dlt = trapz(f(2:end), H(2:end)./f(2:end));
    end
    y_W = 1/C_dlt * real(trapz(s, W./s.^1.5)*2);
    y_T = 1/C_dlt * real(trapz(s, T./s.^1.5)*2);

    %%
    mserror_W = mean((x - y_W).^2);
    mserror_T = mean((x - y_T).^2);

    figure(1);
    subplot(2, 1, 1);
    plot(x);
    hold on
    plot(y_W, '-.m');
    hold off
    title(['K = ', num2str(K, '%d'), ', Mean suqare error: ', num2str(mserror_W, '%e')]);
    subplot(2, 1, 2);
    plot(x);
    hold on
    plot(y_W, '-.m');
    hold off
    xlim([2000, 2200]);

    figure(2);
    subplot(2, 1, 1);
    plot(x);
    hold on
    plot(y_T, '-.m');
    hold off
    title(['K = ', num2str(K, '%d'), ', Mean suqare error: ', num2str(mserror_T, '%e')]);
    subplot(2, 1, 2);
    plot(x);
    hold on
    plot(y_T, '-.m');
    hold off
    xlim([2000, 2200]);

    figure(3);
    subplot(2, 2, 1);
    p = pcolor(1:size(W, 2), F, abs(W).^2);
    p.EdgeColor = 'none';
    titlestr = sprintf('K = %d', K);
    for i=1:numel(f0)
        titlestr = [titlestr, ', f', num2str(i - 1, '%d'), ' = ', num2str(f0(i), '%3.4f')];
    end
    title(titlestr);
    subplot(2, 2, 2);
    p = pcolor(1:size(T, 2), F, abs(T).^2);
    p.EdgeColor = 'none';
    subplot(2, 2, 3);
    p = pcolor(1:size(W, 2), F, log(abs(W).^2));
    p.EdgeColor = 'none';
    subplot(2, 2, 4);
    p = pcolor(1:size(T, 2), F, log(abs(T).^2));
    p.EdgeColor = 'none';

    figure(4);
    plot(f, morsewavelet(gam, be, 0, f));
    hold on
    plot(f, H);
    hold off
end