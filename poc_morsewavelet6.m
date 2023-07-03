function poc_morsewavelet6
    % data simulation
    fs = 2000;
    f0 = (fs*0.3)*rand;
    t = (1:1500)./fs;
    a = sqrt(1./(1:1500));
    x = a.*sin(2*pi*f0.*t);
    
    addpath(strcat(userpath, '/lib2/numdiff/'));

    % test
    be = 30;
    gam = 9;
    f = linspace(0, 1, numel(x));
    f_c = (be/gam)^(1/gam)/(2*pi);
    s = f_c/(f0/fs);
    X = fft(x);

    r = rand(2, 1);
    r = r./sum(r);
    [H_0, dH_0] = morsewavelet(gam, be, 0, s.*f);
    [H_1, dH_1] = morsewavelet(gam, be, 1, s.*f);

    W = ifft(X.*(r(1).*H_0 + r(2).*H_1));
    dW = ifft(X.*(r(1).*dH_0 + r(2).*dH_1)./s);

    Omg = -1i.*dW./(2*pi.*W);

    % plot
    figure(1);
    subplot(4, 1, 1);
    plot(t, abs(W));
    hold on
    plot(t, a);
    hold off
    title(['f_0 = ', num2str(f0, '%3.4f')]);

    subplot(4, 1, 2);
    D = numdiff.dif(real(W), 1);
    plot(t, D);
    hold on
    plot(t, real(dW), '-.m');
    hold off

    subplot(4, 1, 3);
    D = numdiff.dif(imag(W), 1);
    plot(t, D);
    hold on
    plot(t, imag(dW), '-.m');
    hold off

    subplot(4, 1, 4);
    plot(t, fs.*real(Omg));
    hold on;
    plot([t(1), t(end)], [f0, f0], '-.m');
    hold off
end