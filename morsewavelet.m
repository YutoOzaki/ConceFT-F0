function [H, dH, ddH] = morsewavelet(gam, be, k, f)
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
    r = (2*be + 1)/gam;
    c = r - 1;
    A = sqrt(pi*gam*2^r*gamma(k + 1)/gamma(k + r));
    x = 2*(2*pi.*f).^gam;
    m = (0:k)';
    L = sum(x.^m .* (-1).^m.*gamma(k + c + 1)./(gamma(c + m + 1).*gamma(k - m + 1))./factorial(m), 1);
    H = sqrt(2)*A*(2*pi*f).^be.*exp(-(2*pi.*f).^gam).*L;

    %% Derivative for instantenous frequency estimation
    dH = H.*f;
    ddH = dH.*f;
end