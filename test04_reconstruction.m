function test04_reconstruction
    %%
    fs = 400 + round(200*rand);
    N = 4000 + round(500*rand);
    
    %%{
    sgm = 3*rand;
    x = normrnd(0, sgm, [1, N]);
    %}
    
    %{
    sgm = 2*rand;
    x = (sqrt(2)*sgm).*sin((2 + 60*rand)*2*pi.*(0:N - 1)./fs);
    %}
    
    %%
    be = 20 + rand*10;
    gam = 1 + rand*9;
    k = 0;
    f = linspace(0, 1, numel(x));
    X = fft(x);
    [~, ome_c] = morsefreq(be, gam);
    [C, C_dlt] = morseadmissibility(gam, be);

    %%
    s_min = ome_c/0.5;
    s_max = ome_c/(1/numel(x));
    J = log2(s_max) - log2(s_min);
    dj = 0.01 + rand*0.15;
    s = s_min*2.^(0:dj:J)';
    
    %%
    W = zeros(numel(s), numel(x));
    E = zeros(numel(s), 1);
    for i=1:numel(s)
        H = morsewavelet(gam, be, k, s(i).*f);
        W(i, :) = ifft(X.*H);
        E(i) = trapz(s(i).*f, H.^2);
    end
    W = sqrt(s).*W;

    %% Reconstruction using a delta function as the synthesis wavelet
    y = 1/C_dlt * real(trapz(s, W./s.^1.5)*2);

    %% Energy conservation (Parseval's theorem)
    E_W = 1/C * trapz(s, sum(abs(W).^2./s.^2, 2)/N)*2;
    E_x = sum(x.^2)/numel(x);
    fprintf('E_x/E_W = %3.6f\n', E_x/E_W);
    fprintf('sgm^2/E_W = %3.6f\n', sgm^2/E_W);

    %% plot
    F = ome_c./s.*fs;
    t = (0:(numel(x) - 1))./fs;

    figure(1);
    subplot(2, 1, 1);
    p = pcolor(t, s, abs(W).^2);
    p.EdgeColor = 'none';
    subplot(2, 1, 2);
    p = pcolor(t, F, abs(W).^2);
    p.EdgeColor = 'none';

    sqerror = mean((x - y).^2);
    figure(2);
    subplot(2, 1, 1);
    plot(x); hold on;
    plot(y, '-.m'); hold off
    title(['Squared error = ', num2str(sqerror, '%e')]);
    subplot(2, 1, 2);
    plot(x); hold on;
    plot(y, '-.m'); hold off
    xlim([2000, 2100]);
end