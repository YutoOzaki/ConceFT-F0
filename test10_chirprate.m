function test10_chirprate
    %% data simulation
    N = 1600;
    fs = 450;
    ph_0 = rand*2*pi;
    f0 = 5 + 10*rand;
    f1 = 70 + 20*rand;
    t = (0:(N - 1))./fs;

    %{
    c = (f1 - f0)/(2 + 2*rand);
    f_x = f0 + t.*c;
    h_f = @(i) f0 + c*(i - 1)/fs;
    h_q = @(t) c.*ones(numel(t), 1);
    x = cos(ph_0 + 2*pi.*(c/2.*t.^2 + f0.*t));
    %}

    %{
    c = (f1 - f0)/(2 + 2*rand);
    f_x = f0 + 0.5.*t.^2.*c;
    h_f = @(i) f0 + 0.5*c*((i - 1)/fs)^2;
    h_q = @(t) c.*t;
    x = sin(ph_0 + 2*pi.*(0.5*(1/3)*c.*t.^3 + f0.*t));
    %}
    
    %%{
    c = (f1/f0)^(1/(N/fs));
    f_x = f0.*c.^t;
    h_f = @(i) f0*c^((i - 1)/fs);
    h_q = @(t) t.*f0.*c.^(t - 1);
    x = sin(ph_0 + 2*pi*f0.*(c.^t - 1)./log(c));
    %}

    figure(1);
    clf; cla;
    spectrogram(x, hann(32), 31, 32, fs, 'yaxis');
    colorbar off
    hold on; plot(t, f_x, 'Color', 'r'); hold off;

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

    X = fft(x);
    n = (0:numel(x) - 1);%./numel(x);
    Omg = zeros(numel(s), numel(x));
    tau = zeros(numel(s), numel(x));
    q = zeros(numel(s), numel(x));
    qhat = zeros(numel(s), numel(x), 2);
    fhat = zeros(numel(x), 1);
    chat = zeros(numel(x), 3);
    gdhat = zeros(numel(x), 1);
    E = zeros(numel(s), 1);
    for i=1:numel(s)
        [H, xiH, xisqH, dH, xidH] = morsewavelet(gam, be, k, s(i).*f);
        E(i) = trapz(s(i).*f(1:N/2), H(1:N/2).^2);

        W_H = ifft(X.*H).*sqrt(s(i));
        W_xiH = ifft(X.*xiH).*sqrt(s(i));
        W_xisqH = ifft(X.*xisqH).*sqrt(s(i));
        W_dH = ifft(X.*dH).*sqrt(s(i));
        W_xidH = ifft(X.*xidH).*sqrt(s(i));

        Omg(i, :) = (1/s(i)).*W_xiH./W_H;
        tau(i, :) = n + s(i)/(1i*2*pi).*(W_dH./W_H);
        q(i, :) = (1i*2*pi)/s(i)^2 .* (W_H.*W_xisqH - W_xiH.^2) ./ (W_H.^2 + W_H.*W_xidH - W_dH.*W_xiH);
        qhat(i, :, 1) = real(gradient(Omg(i, :), 1)./gradient(tau(i, :), 1));
    end
    q = real(q);

    for i=1:numel(x)
        qhat(:, i, 2) = real(gradient(Omg(:, i), 1)./gradient(tau(:, i), 1));

        f_i = h_f(i);
        s_i = (f_i/fs/omg_c)^-1;
        [~, idx] = min(abs(s_i - s));

        fhat(i) = Omg(idx, i);
        gdhat(i) = tau(idx, i);
        chat(i, 1) = q(idx, i);
        chat(i, 2) = qhat(idx, i, 2);
        chat(i, 3) = qhat(idx, i, 1);
    end
    fhat = real(fhat);
    gdhat = real(gdhat);
    chat = real(chat);

    %% validation
    figure(2);
    subplot(2, 1, 1); plot(E);
    subplot(2, 1, 2); plot(abs(E - 1) < 1e-8);
    
    figure(3);
    plot(t, fhat.*fs);
    hold on
    plot(t, f_x, '-.m');
    hold off
    title(sprintf('error = %e', mean(abs(f_x(:) - fhat(:).*fs))));
    
    figure(4);
    plot(t, gdhat)
    hold on
    plot(t, 1:N, '-.m');
    hold off
    title(sprintf('error = %e', mean(abs(t(:) - gdhat(:).*(N/fs)))));

    figure(5);
    subplot(2, 1, 1);
    plot(t, chat(:, 1).*fs^2);
    hold on
    plot(t, chat(:, 2).*fs^2, 'Color', "#D95319");
    plot(t, chat(:, 3).*fs^2, 'Color', "#EDB120", 'LineStyle', '--');
    hold off
    subplot(2, 1, 2);
    plot(t, h_q(t), '-.m');
    yl = ylim();
    hold on
    plot(t, chat(:, 2).*fs^2, 'Color', "#D95319");
    plot(t, chat(:, 3).*fs^2, 'Color', "#EDB120", 'LineStyle', '--');
    hold off
    ylim(yl);
end