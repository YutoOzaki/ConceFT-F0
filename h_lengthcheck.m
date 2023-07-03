function L_integer = h_lengthcheck(be, gam, F, fs)
    %%
    f_c = (be/gam)^(1/gam)/(2*pi);
    s = f_c/(F/fs);

    %%
    L = fzero(@(N) h_unitenergycheck(gam, be, s, round(N)), [round(fs/F*0.5), fs*2]);
    L_integer = round(L);
end

function d = h_unitenergycheck(gam, be, s, N)
    %%
    eps = 1e-4;

    %%
    f = s.*linspace(0, 1, N);
    H = morsewavelet(gam, be, 0, f);
    E = trapz(f, H.^2);

    %%
    d = abs(1 - E) - eps;
end