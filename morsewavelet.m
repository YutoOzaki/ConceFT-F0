function [H, dH, ddH, gH, gdH] = morsewavelet(gam, be, k, f)
    %% Reference
    % [1] Olhede, S., C & Walden, A., T. (2002). Generalized Morse Wavelets.
    % IEEE Transactions on Signal Processing, 50(11), 2661-2670.

    %% Memo
    % The support of f is [0, 1] (normalized frequency), like
    % linspace(0, 1, 256). Also, see that f is multiplied by 2π in the
    % implementation below.
    %
    % The center frequency of the 0th order Morse wavelet is
    % (be/gam)^(1/gam)/(2*pi). Again, this value lies in the range of [0, 1].
    %
    % This Morse wavelet has unit energy, which can be confirmed by 
    % trapz(f, H.^2) being equal to 1.

    %% Equation (10) from [1].
    xi = (2*pi).*f;

    r = (2*be + 1)/gam;
    c = r - 1;
    A = sqrt(pi*gam*2^r*gamma(k + 1)/gamma(k + r));
    P = xi.^be;
    Q = exp(-xi.^gam);
    L = laguerrepoly(2.*xi.^gam, k, c);
    H = sqrt(2)*A*P.*Q.*L;

    %% Derivative for instantenous frequency estimation
    dH = H.*f;
    ddH = dH.*f;
    
    %% Differentiation
    dP = be*(2*pi)^be.*f.^(be - 1);
    dQ = -gam*(2*pi)^gam.*f.^(gam - 1).*Q;
    dL = 2*gam*(2*pi)^gam.*f.^(gam - 1).*(-laguerrepoly(2.*xi.^gam, k - 1, c + 1));
    
    gH = sqrt(2)*A.*(...
        dP.*Q.*L + ...
        P.*dQ.*L + ...
        P.*Q.*dL);
    
    gdH = (gH.*f + H);
end