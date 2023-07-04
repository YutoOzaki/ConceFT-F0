function test08_gradient
    %%
    addpath('C:\Users\yuto\Documents\MATLAB\lib2\numdiff');
    N = 512;
    x = rand(1, N);
    X = fft(x);

    %%
    k = randi(5) - 1;
    be = 20;
    gam = 12;

    f = linspace(0, 1, N);
    s = rand*2;

    [H, dH, ddH] = morsewavelet(gam, be, k, s.*f);
    E = trapz(s.*f, H.^2);
    E_d = trapz(s.*f, dH.^2);
    W = ifft(X.*H).*sqrt(s);
    WW = ifft(X.*dH).*sqrt(s);

    %% analytical
    dW = (2*1i*pi/sqrt(s)).*ifft(X.*dH);
    dWW = (2*1i*pi/sqrt(s)).*ifft(X.*ddH);
    Omg = (1/s).*WW./W;

    %% numerical
    A = numdiff.dif(W, 1);
    B = numdiff.dif(WW, 1);
    C = (1/(2*1i*pi)).*numdiff.dif(W, 1)./W;
    
    %% error
    re_r = abs(real(dW) - real(A))./real(A);
    re_i = abs(imag(dW) - imag(A))./imag(A);

    re_r_B = abs(real(dWW) - real(B))./real(B);
    re_i_B = abs(imag(dWW) - imag(B))./imag(B);
    
    re_O = abs(imag(Omg) - imag(C))./imag(C);

    %% plot
    figure(1);
    clf; cla;
    subplot(3, 1, 1);
    plot(f, H);
    title(sprintf('k = %d, s = %3.4f, E = %e', k, s, E));
    subplot(3, 1, 2);
    plot(real(dW)); hold on; plot(real(A), '-.m'); hold off
    title(sprintf('E[re_r] = %e', mean(re_r)));
    subplot(3, 1, 3);
    plot(imag(dW)); hold on; plot(imag(A), '-.m'); hold off
    title(sprintf('E[re_i] = %e', mean(re_i)));

    figure(2);
    clf; cla;
    subplot(3, 1, 1);
    plot(f, dH);
    title(sprintf('k = %d, s = %3.4f, E = %e', k, s, E_d));
    subplot(3, 1, 2);
    plot(real(dWW)); hold on; plot(real(B), '-.m'); hold off
    title(sprintf('E[re_r] = %e', mean(re_r_B)));
    subplot(3, 1, 3);
    plot(imag(dWW)); hold on; plot(imag(B), '-.m'); hold off
    title(sprintf('E[re_i] = %e', mean(re_i_B)));

    figure(3);
    plot(real(Omg)); hold on; plot(real(C), '-.m'); hold off
    title(sprintf('E[re_O] = %e', mean(re_O)));
end