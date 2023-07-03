function test03_gaussiannoise
    %% data generation
    sgm = rand*2;
    x = normrnd(0, sgm, [1, 4096]);

    %% test
    be = 30;
    gam = 9;
    K = 4;
    f = linspace(0, 1, numel(x));
    X = fft(x);

    s = rand(2048, 1).*3;
    sgmsq_W = zeros(numel(s), K + 1);
    E = zeros(numel(s), K + 1);

    for k=0:K
        for i=1:numel(s)
            H = morsewavelet(gam, be, k, s(i).*f);
            E(i, k + 1) = trapz(s(i).*f, H.^2);
            W = sqrt(s(i)).*ifft(X.*H);
            sgmsq_W(i, k + 1) = mean(abs(W).^2);
        end
    end
    
    %% plot
    figure(1);
    for k=0:K
        subplot(K + 1, 1, k + 1);
        idx = abs(E(:, k + 1) - 1) < 1e-8;
        histogram(sqrt(sgmsq_W(idx, k + 1)));
        sgmsq_W_k = mean(sqrt(sgmsq_W(idx, k + 1)));
        hold on
        yl = ylim();
        plot([sgm, sgm], yl, 'Color', 'k');
        plot([sgmsq_W_k, sgmsq_W_k], yl, 'Color', 'm');
        hold off
        title(['k = ', num2str(k, '%d')]);
    end
    
    figure(2);
    for k=0:K
        subplot(K + 1, 1, k + 1);
        idx = abs(E(:, k + 1) - 1) < 1e-8;
        scatter(s(idx), sqrt(sgmsq_W(idx, k + 1)), 'Marker', '.');
        hold on
        plot([min(s(idx)), max(s(idx))], [sgm, sgm], 'Color', 'k');
        hold off
        title(['k = ', num2str(k, '%d')]);
    end

    figure(3);
    subplot(2, 1, 1);
    scatter(s, mean(sqrt(sgmsq_W), 2), 'Marker', '.');
    hold on
    plot([min(s), max(s)], [sgm, sgm], 'Color', 'k');
    hold off
    yl = ylim();

    subplot(2, 1, 2);
    idx = 0;
    for k=0:K
        idx = idx | abs(E(:, k + 1) - 1) < 1e-8;
    end
    scatter(s(idx), mean(sqrt(sgmsq_W((idx), :)), 2), 'Marker', '.');
    hold on
    plot([min(s(idx)), max(s(idx))], [sgm, sgm], 'Color', 'k');
    hold off
    ylim(yl);
end