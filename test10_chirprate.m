function test10_chirprate
    %% data simulation
    N = 512;
    fs = 200;
    f0 = 5 + 10*rand;
    f1 = 80 + 10*rand;
    t = (0:(N - 1))./fs;
    ph_0 = rand*2*pi;
    c = (f1 - f0)/t(end);
    x = sin(ph_0 + 2*pi.*(c/2.*t.^2 + f0.*t));

    figure(1);
    clf; cla;
    spectrogram(x, hann(32), 31, 32, fs, 'yaxis');
    colorbar off
    hold on; plot([t(1), t(end)], [f0, f1], 'Color', 'r'); hold off;

    %% CWT and chirp rate
    gam = 9;
    be = 17;
    k = 0;
    f = linspace(0, 1, numel(x));

    [~, omg_c] = morsefreq(be, gam);
    s_min = omg_c/0.5;
    s_max = omg_c/(1/numel(x));
    J = log2(s_max) - log2(s_min);
    dj = 1/32;
    s = s_min*2.^(0:dj:J)';

    addpath('C:\Users\yuto\Documents\MATLAB\lib2\numdiff');
    X = fft(x);
    Omg = zeros(numel(s), numel(x));
    tau = zeros(numel(s), numel(x));
    q = zeros(numel(s), numel(x));
    qhat = zeros(numel(s), numel(x));
    E = zeros(numel(s), 1);
    for i=1:numel(s)
        [H, xiH, xisqH] = morsewavelet(gam, be, k, s(i).*f);
        E(i) = trapz(s(i).*f(1:N/2), H(1:N/2).^2);
        
        if abs(E(i) - 1) < 1e-8 
            dH = numdiff.dif(H, 1);
            xidH = numdiff.dif(xiH, 1);
    
            W_H = ifft(X.*H).*sqrt(s(i));
            W_xiH = ifft(X.*xiH).*sqrt(s(i));
            W_xisqH = ifft(X.*xisqH).*sqrt(s(i));
            W_dH = ifft(X.*dH).*sqrt(s(i));
            W_xidH = ifft(X.*xidH).*sqrt(s(i));
    
            Omg(i, :) = (1/s(i)).*W_xiH./W_H;
            tau(i, :) = t + s(i)/(1i*2*pi).*(W_dH./W_H);
            q(i, :) = (1i*2*pi)/s(i)^2 .* (W_H.*W_xisqH - W_xiH.^2) ./ (W_H.^2 + W_H.*W_xidH + W_dH.*W_xiH);
        end
    end
    q = real(q);
    
    for i=1:numel(x)
        qhat(:, i) = real(numdiff.dif(Omg(:, i), 1)./numdiff.dif(tau(:, i), 1));
    end
    qhat = real(qhat);

    %% validation
    figure(2);
    subplot(2, 1, 1); plot(E);
    subplot(2, 1, 2); plot(abs(E - 1) < 1e-8);
    
    i = randi(100) + 10;
    f_s = omg_c/s(i)*fs;
    t_s = (f_s - f0)/c;
    [~, idx] = min(abs(t - t_s));
    
    figure(3);
    clf; cla;
    plot(t, fs.*q(i, :));
    hold on;
    plot(t, fs.*qhat(i, :));
    yl = ylim(); plot([t_s, t_s], yl, '-.m');
    hold off
    title(sprintf('f = %3.4f, IF = %3.4f, c = %3.4f, qhat = %3.4f, q = %3.4f',...
        f0 + c*t_s, fs*real(Omg(i, idx)), c, fs*qhat(i, idx), fs*q(i, idx)));

    figure(4);
    plot(t, real(tau(i, :) - t));
    hold on;
    yl = ylim(); plot([t_s, t_s], yl, '-.m');
    hold off
    title(sprintf('dt = %3.4f, dtau = %3.4f', t(idx) - t_s, real(tau(i, idx) - t(idx))));
end