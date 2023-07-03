function [C, C_dlt, a] = morseadmissibility(gam, be)
    %% Reference related to C_dlt
    % [1] Farge, M. (1992). Wavelet transforms and their applications to
    % turbulennce. Annu. Rev. Fluid Mech., 24, 395-457.
    % [2] Torrence, C. & Compo, G. P. (1998). A Practical Guide to Wavelet
    % Analysis. Bulletin of the American Meteorological Society, 79(1),
    % 61-78.
    % [3] Holschneider, M. & Tchamitchian, Ph. (1991). Pointwise analysis
    % of Riemann's "nondifferentiable" function. Inventiones mathematicae,
    % 105, 157-175.

    %% Reference related to C (admissibility constant)
    % [4] Lilly, J. M. & Olhede, S. C. (2009). Higher-Order Properties of
    % Analytic Wavelets. IEEE Transactions on Signal Processing, 57(1),
    % 146-160.

    %% Memo
    % These constants are for the case of 0th order Morse wavelets (k = 0).
    % C_dlt is used for reconstruction using delta function as a synthesis
    % wavelet.

    %%
    k = 0;

    r = (2*be + 1)/gam;
    A = sqrt(pi*gam*2^r*gamma(k + 1)/gamma(k + r));
    a = sqrt(2)*A;

    %%
    C = a^2/(gam*2^(2*be/gam)) * gamma(2*be/gam);
    C_dlt = a/gam * gamma(be/gam);
end