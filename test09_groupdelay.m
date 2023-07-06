function test09_groupdelay
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

    %% CWT and group delay
    gam = 9;
    be = 17;
    k = 0;
    f = linspace(0, 1, numel(x));
    s = rand*1.2;

    X = fft(x);
    H = morsewavelet(gam, be, k, s.*f);
    W_H = ifft(X.*H).*sqrt(s);
    addpath('C:\Users\yuto\Documents\MATLAB\lib2\numdiff');
    dH = numdiff.dif(H, 1);
    W_dH = ifft(X.*dH).*sqrt(s);
    
    dtau = s/(1i*2*pi).*W_dH./W_H;
    tau = t + dtau;
    
    %% validation
    [~, omg_c] = morsefreq(be, gam);
    f_s = omg_c/s*fs;
    t_s = (f_s - f0)/c;
    [~, idx] = min(abs(t - t_s));

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