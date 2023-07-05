function [H, dH, ddH, gH] = morsewavelet(gam, be, k, f)
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
    x = 2.*xi.^gam;
    m = (0:k)';
    C_m = (-1).^m.*gamma(k + c + 1)./(gamma(c + m + 1).*gamma(k - m + 1))./factorial(m);
    L = sum(C_m .* x.^m, 1);
    H = sqrt(2)*A*xi.^be.*exp(-xi.^gam).*L;

    %% Derivative for instantenous frequency estimation
    dH = H.*f;
    ddH = dH.*f;
    
    %% Differentiation
    df = f(end)/(numel(f) - 1);
    %gL = sum(C_m .* gam.*m.*pi.^(gam.*m).*2.^((gam + 1).*m).*(x.^gam).^m./x, 1);
    m = (0:k - 1)';
    c = r;
    C_m = (-1).^m.*gamma(k + c + 1)./(gamma(c + m + 1).*gamma(k - m + 1))./factorial(m);
    gL = -sum(C_m .* x.^m, 1) .* (gam*pi^gam*2^(gam + 1)*f.^(gam - 1));

    gH = sqrt(2)*A * ((be*xi.^(be - 1)).*exp(-xi.^gam).*L + xi.^be.*(-gam.*xi.^(gam - 1).*exp(-xi.^gam)).*L + ...
        xi.^be.*exp(-xi.^gam).*gL).*(2*pi*df);
    gH(1) = 0;
end