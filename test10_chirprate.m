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
    q = zeros(numel(s), numel(x));
    E = zeros(numel(s), 1);
    for i=1:numel(s)
        [H, xiH, xisqH] = morsewavelet(gam, be, k, s(i).*f);
        E(i) = trapz(s(i).*f(1:N/2), H(1:N/2).^2);
        
        if abs(E(i) - 1) < 1e-8 
            dH = numdiff.dif(H, s(i));
            xidH = numdiff.dif(xiH, s(i));
    
            W_H = ifft(X.*H).*sqrt(s(i));
            W_xiH = ifft(X.*xiH).*sqrt(s(i));
            W_xisqH = ifft(X.*xisqH).*sqrt(s(i));
            W_dH = ifft(X.*dH).*sqrt(s(i));
            W_xidH = ifft(X.*xidH).*sqrt(s(i));
    
            q(i, :) = (1i*2*pi)/s(i)^2 .* (W_H.*W_xisqH - W_xiH.^2) ./ (W_H.^2 + W_H.*W_xidH + W_dH.*W_xiH);
        end
    end
    q = real(q);

    %% validation
    figure(2);
    subplot(2, 1, 1); plot(E);
    subplot(2, 1, 2); plot(abs(E - 1) < 1e-8);
    
    i = 40;
    f_s = omg_c/s(i)*fs;
    t_s = (f_s - f0)/c;
    [~, idx] = min(abs(t - t_s));

    [H, xiH, ~] = morsewavelet(gam, be, k, s(i).*f);
    dH = numdiff.dif(H, s(i));
    W_H = ifft(X.*H).*sqrt(s(i));
    W_xiH = ifft(X.*xiH).*sqrt(s(i));
    W_dH = ifft(X.*dH).*sqrt(s(i));

    Omg = (1/s(i)).*W_xiH./W_H;
    tau = t + s(i)/(1i*2*pi).*(W_dH./W_H);
    qhat = real(numdiff.dif(Omg, 1)./numdiff.dif(tau, 1));
    
    figure(3);
    clf; cla;
    plot(t, q(i, :));
    hold on;
    plot(t, qhat);
    yl = ylim(); plot([t_s, t_s], yl, '-.m');
    hold off

    figure(2);
    clf; cla;
    plot(t, real(tau));
    hold on;
    yl = ylim(); plot([t_s, t_s], yl, '-.m');
    xl = xlim(); plot(xl, real([tau(idx), tau(idx)]), '-.b');
    hold off;
    
    figure(3);
    subplot(2, 1, 1);
    plot(f, H);
    E = trapz(f(1:N/2).*s, H(1:N/2).^2);
    title(sprintf('E = %e', E));
    subplot(2, 1, 2);
    plot(f, dH);

    figure(4);
    clf; cla;
    plot(t, real(dtau));
    hold on;
    yl = ylim(); plot([t_s, t_s], yl, '-.m');
    xl = xlim(); plot(xl, real([dtau(idx), dtau(idx)]), '-.b');
    hold off;
    
    figure(5);
    subplot(2, 1, 1);
    plot(t, abs(W_H).^2);
    hold on; yl = ylim(); plot([t_s, t_s], yl, '-.m'); hold off
    title(sprintf('f_W = %3.4f [Hz]', f_s));
    subplot(2, 1, 2);
    plot(t, abs(W_dH).^2);
    hold on; yl = ylim(); plot([t_s, t_s], yl, '-.m'); hold off
end