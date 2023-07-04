function poc_morsewavelet2
    N = 512;

    be = 30;
    gam = 9;
    f = linspace(0, 1, N);
    k = 0;
    H = morsewavelet(gam, be, k, f);
    E = trapz(f, H.^2);
    omg_c = (be/gam)^(1/gam)/(2*pi);

    s = 1.37;
    H_s = sqrt(s).*morsewavelet(gam, be, k, s.*f);
    E_s = trapz(f, H.^2);
    omg_s = omg_c/s;
    
    figure(1);
    subplot(2, 1, 1);
    plot(f, H);
    hold on
    yl = ylim;
    plot(omg_c.*[1, 1], yl);
    hold off
    title(['Energy = ', num2str(E, '%3.4f')]);

    subplot(2, 1, 2);
    plot(f, H_s);
    hold on
    yl = ylim;
    plot(omg_s.*[1, 1], yl);
    hold off
    title(['Energy = ', num2str(E_s, '%3.4f')]);
end