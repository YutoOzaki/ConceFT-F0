function test09_groupdelay
    rng(11);

    %%
    addpath('C:\Users\yuto\Documents\MATLAB\lib2\numdiff');
    N = 512;
    x = rand(1, N);
    X = fft(x);

    %%
    %k = randi(5) - 1;
    k = 1;
    be = 20;
    gam = 12;

    f = linspace(0, 1, N);
    s = rand*2;

    [H, ~, ~, kH] = morsewavelet(gam, be, k, s.*f);

    %% numerical
    A = numdiff.dif(H, 1);
    re = abs(kH - A)./A;
    
    %%
    E = trapz(s.*f, H.^2);
    E_k = trapz(s.*f, kH.^2);
    W = ifft(X.*H).*sqrt(s);
    kW = ifft(X.*kH).*sqrt(s);
    
    dtau = s/(2*1i*pi).*kW./W;

    %%
    figure(1);
    clf; cla;
    subplot(4, 1, 1);
    plot(f, H);
    title(sprintf('k = %d, s = %3.4f, E = %e', k, s, E));
    subplot(4, 1, 2);
    plot(f, kH);
    title(sprintf('k = %d, s = %3.4f, E = %e', k + 1, s, E_k));
    subplot(4, 1, 3);
    plot(real(dtau)); hold on; plot(imag(dtau), 'Linestyle', '-.'); hold off
    subplot(4, 1, 4);
    plot(f, kH); hold on; plot(f, A, '-.m'); hold off
    title(sprintf('re = %e', mean(re, "omitnan")));
end