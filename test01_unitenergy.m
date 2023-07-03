function test01_unitenergy
    %%
    N = 64000;
    be = 3 + rand(32, 1).*20;
    gam = 1 + rand(32, 1).*(2.*be + 1);

    %%
    K = 7;
    f = linspace(0, 1, N);
    E = zeros(numel(be), K + 1);

    for k=0:K
        for i=1:numel(be)
            H = morsewavelet(gam(i), be(i), k, f);
            E(i, k + 1) = trapz(f, H.^2);
        end
    end
    
    assert(all(abs(E(:) - 1) < 1e-8), 'Check unit energy normalization');
    
    %%
    idx = randperm(numel(be), K + 1);

    figure;
    for k=0:K
        H = morsewavelet(gam(idx(k + 1)), be(idx(k + 1)), k, f);
        subplot(4, 2, k + 1);
        plot(f, H);
        title(['k = ', num2str(k, '%d'), ', \beta = ', num2str(be(idx(k + 1)), '%3.3f'),...
            ', \gamma = ', num2str(gam(idx(k + 1)), '%3.3f'), ', E = ', num2str(E(idx(k + 1), k + 1), '%e')]);
    end
end