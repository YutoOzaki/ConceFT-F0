function U = h_gausscorr(P, voice)
    %%
    M = size(P, 1);
    I = 1:M;
    U = P .* 0;
    J = 8;
    %r = 1./I;
    
    %%
    for m=1:M
        idx_u = m + voice*log2(1:J);
        idx_d = [m - voice, (idx_u(2:end) + idx_u(1:end - 1))./2];
        idx_0 = [(idx_u + idx_d)./2, (idx_u(1:J - 1) + idx_d(2:end))./2];

        K = zeros(1, M);
        theta_ini = pi;
        for i=1:(J - 1)
            idx_p = find(idx_d(i) <= I, 1, 'first');
            idx_q = find(I <= idx_d(i + 1), 1, 'last');
            theta = spline([idx_d(i), idx_0(i), idx_u(i), idx_0(J + i), idx_d(i + 1)],...
                theta_ini + [0, 0.5*pi, pi, 1.5*pi, 2*pi], idx_p:idx_q);
            K(idx_p:idx_q) = cos(theta);
        end
        idx_q = find(I >= idx_0(end), 1, 'first');
        K(idx_q:end) = 0;
        idx_p = find(I < idx_d(1), 1, 'last');
        K(1:idx_p) = -1;

        %K = K.*r;

        U(m, :) = K*P;
    end
end