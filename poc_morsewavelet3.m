function poc_morsewavelet3
    fs = 2000;
    a = rand(3, 1);
    f0 = sort((fs*0.3).*rand(3, 1));
    t = (1:1500)./fs;
    x = 0;
    for i=1:numel(a)
        x = x + a(i).*sin(2*pi*f0(i).*t);
    end
    
    be = 30;
    gam = 9;
    k = 0;
    omg_c = (be/gam)^(1/gam)/(2*pi);
    f = linspace(0, 1, numel(x));
    W = zeros(numel(a), numel(x));
    H = zeros(numel(a), numel(x));
    E = zeros(numel(a), 1);
    X = fft(x);

    for i=1:size(W, 1)
        s = omg_c/(f0(i)/fs);
        W(i, :) = ifft(X.*morsewavelet(gam, be, k, s.*f));
        H(i, :) = morsewavelet(gam, be, k, s.*f)./morsewavelet(gam, be, k, omg_c).*2;
        E(i) = trapz(f, H(i, :).^2);
    end
    W = W./morsewavelet(gam, be, k, omg_c).*2;
    re = abs(a - abs(W(:, 500)))./a;

    figure(1);
    subplot(2, 1, 1);
    titlestr = 'Relative error = ';
    for i=1:numel(a)
        plot(t, abs(W(i, :)), 'Color', 'b');
        hold on
        plot([t(1), t(end)], [a(i), a(i)], '-.m');
        titlestr = [titlestr, num2str(re(i), '%e'), ', '];
    end
    hold off
    xlim([t(1), t(end)]);
    title(titlestr);

    subplot(2, 1, 2);
    titlestr = 'Energy = ';
    for i=1:numel(a)
        plot(f, H(i, :), 'Color', 'b');
        hold on
    end
    yl = ylim;
    for i=1:numel(a)
        plot((f0(i)/fs).*[1, 1], yl, '-.m');
        titlestr = [titlestr, num2str(E(i), '%3.4f'), ', '];
    end
    hold off
    title(titlestr)
end