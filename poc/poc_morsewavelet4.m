function poc_morsewavelet4
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

    %{
    omg = s.*(2*pi).*f;
    al = 2*(exp(1)*gam/be)^(be/gam);
    H = al*omg.^be.*exp(-omg.^gam);
    W = ifft(X.*H);
    
    al_b = 2*(exp(1)*gam/(be + 1))^((be + 1)/gam);
    C = (-1/1i)*al/al_b./s;
    dH = C .* al_b*omg.^(be + 1).*exp(-omg.^gam);
    dW = ifft(X.*dH);
    %}

    %%{
    k = 1;
    [H, dH] = morsewavelet(gam, be, k, s.*f);
    %C = 1/morsewavelet(gam, be, k, f_c).*2;
    C = 1;
    W = ifft(X.*H).*C;
    dH = dH./s;
    dW = ifft(X.*dH).*C;
    %}

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