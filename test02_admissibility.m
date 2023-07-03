function test02_admissibility
    %%
    N = 64000;
    be = 3 + rand(32, 1).*20;
    gam = 1 + rand(32, 1).*(2.*be + 1);
    k = 0;

    %%
    f = linspace(0, 1, N);
    C = zeros(numel(be), 3);
    C_dlt = zeros(numel(be), 2);
    
    for i=1:numel(be)
        H = morsewavelet(gam(i), be(i), k, f);
        E = trapz(f, H.^2);
        
        if abs(E - 1) < 1e-8
            [C_i, C_dlt_i] = morseadmissibility(gam(i), be(i));

            C(i, 1) = C_i;
            C(i, 2) = trapz((2*pi).*f(2:end), H(2:end).^2./(2*pi.*f(2:end)));
            C(i, 3) = trapz(f(2:end), H(2:end).^2./f(2:end));

            C_dlt(i, 1) = C_dlt_i;
            C_dlt(i, 2) = trapz(f(2:end), H(2:end)./f(2:end));
        end
    end
    
    assert(all(abs(C(:, 1) - C(:, 2)) < 1e-8), 'Check admissibility constant');
    assert(all(abs(C(:, 1) - C(:, 3)) < 1e-8), 'Check admissibility constant');
    assert(all(abs(C_dlt(:, 1) - C_dlt(:, 2)) < 1e-8), 'Check admissibility constant');

    %%
    figure;
    for i=1:4
        H = morsewavelet(gam(i), be(i), k, f);
        subplot(4, 1, i);
        plot(f, H);
        title(['k = ', num2str(k, '%d'), ', \beta = ', num2str(be(i), '%3.3f'),...
            ', \gamma = ', num2str(gam(i), '%3.3f'), ', C = ', num2str(C(i, 1), '%3.3f')]);
    end
end