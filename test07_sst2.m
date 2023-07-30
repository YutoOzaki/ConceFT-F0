function test07_sst2
    %%
    fs = 2000;
    t = (0:(2*fs - 1))./fs;
    
    %{
    a = rand(3, 1);
    f0 = sort((fs*0.3).*rand(3, 1));
    x = 0;
    for i=1:numel(a)
        x = x + a(i).*sin(2*pi*f0(i).*t);
    end
    f_x = [f0(1).*ones(1, numel(t));...
        f0(2).*ones(1, numel(t));...
        f0(3).*ones(1, numel(t));...
        ];
    %}
    
    %%{
    ph_0 = rand*2*pi;
    f0 = 20 + 10*rand;
    f1 = 900 + 70*rand;
    c = (f1 - f0)/(2 + 2*rand);
    f_x = f0 + 0.5.*t.^2.*c;
    x = sin(ph_0 + 2*pi.*(0.5*(1/3)*c.*t.^3 + f0.*t));
    %}
    
    %{
    ph_0 = rand*2*pi;
    f0 = 20 + 10*rand;
    f1 = 900 + 70*rand;
    c = (f1/f0)^(1/(numel(t)/fs));
    f_x = f0.*c.^t;
    x = sin(ph_0 + 2*pi*f0.*(c.^t - 1)./log(c));
    %}
    
    %{
    f0 = [];
    A = 1 + 3.*t.^2 + 4.*(1 - t).^7;
    ph = 240.*t - 2.*exp(-2.*t).*sin(14*pi.*t);
    x = A.*exp(2*pi.*1i.*ph);
    x = real(x);
    f_x = gradient(ph).*fs;
    %}

    %%
    gam = 2 + rand*9;
    be = 20 + rand*10;
    K = randi(5) - 1;
    %K = 0;
    f = linspace(0, 1, numel(x));
    
    r = rand(K + 1, 1);
    r = r./sum(r);

    %%
    [~, ome_c] = morsefreq(be, gam);
    s_min = ome_c/0.5;
    s_max = ome_c/(1/numel(x));
    J = log2(s_max) - log2(s_min);
    dj = 1/64;
    s = s_min*2.^(0:dj:J)';

    %%
    X = fft(x);
    N = numel(x);
    W = zeros(numel(s), N);
    Omg = zeros(numel(s), N);
    n = 0:(numel(x) - 1);

    for i=1:numel(s)
        H = 0;
        xiH = 0;
        xisqH = 0;
        dH = 0;
        xidH = 0;
        
        for k=0:K
            [H_k, xiH_k, xisqH_k, dH_k, xidH_k] = morsewavelet(gam, be, k, s(i).*f);

            H = H + r(k + 1).*H_k;
            xiH = xiH + r(k + 1).*xiH_k;
            xisqH = xisqH + r(k + 1).*xisqH_k;
            dH = dH + r(k + 1).*dH_k;
            xidH = xidH + r(k + 1).*xidH_k;
        end

        W_H = ifft(X.*H).*sqrt(s(i));
        W_xiH = ifft(X.*xiH).*sqrt(s(i));
        W_xisqH = ifft(X.*xisqH).*sqrt(s(i));
        W_dH = ifft(X.*dH).*sqrt(s(i));
        W_xidH = ifft(X.*xidH).*sqrt(s(i));
        
        W(i, :) = W_H;
        
        Omg_1 = (1/s(i)).*W_xiH./W_H;
        tau = n + s(i)/(1i*2*pi).*(W_dH./W_H);
        
        dOmg = (1i*2*pi/s(i)) .* (W_H.*W_xisqH - W_xiH.^2)./W_H.^2;
        dtau = 1/fs + s(i).*(W_H.*W_xidH - W_dH.*W_xiH)./W_H.^2;
        q = dOmg./dtau;
        Omg_2 = Omg_1 + q.*(n - tau);
        
        Omg(i, :) = Omg_2;
    end
    
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
    
    %% Check energy
    E = 0;
    for j=1:size(f_x, 1)
        for i=1:size(T, 2)
            [~, idx] = min(abs(f_x(j, i) - F));
            E = E + sum(abs(T(idx - 2:idx + 2, i)).^2);
        end
    end
    E = E./sum(abs(T(:)).^2);
    
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
    title(['E = ', num2str(E, '%3.4f')]);
    subplot(2, 2, 3);
    p = pcolor(1:size(W, 2), F, (abs(W).^2).^0.1);
    p.EdgeColor = 'none';
    subplot(2, 2, 4);
    p = pcolor(1:size(T, 2), F, (abs(T).^2).^0.1);
    p.EdgeColor = 'none';

    figure(4);
    plot(f, morsewavelet(gam, be, 0, f));
    hold on
    plot(f, H);
    hold off
end