function [peakfreq, energyfreq] = morsefreq(be, gam)
    peakfreq = (be/gam)^(1/gam) * 1/(2*pi);
    energyfreq = 1/2^(1/gam) * gamma((2*be + 2)/gam)/gamma((2*be + 1)/gam) * 1/(2*pi);
end