function test06_trapz
    M = 1024;
    A = zeros(M, 2);
    c = 1.5;

    for m=1:M
        x = cumsum(rand(1, 32));
        y = rand(1, numel(x));
        A(m, 1) = trapz(x, y./x.^c);
    
        i = randi(numel(x) - 2) + 1;
        j = randi(numel(x) - 2) + 1;
        z = y;
        buf = z(i);
        z(i) = 0;
        f_j = (x(j + 1) - x(j - 1))/x(j)^c;
        f_i = (x(i + 1) - x(i - 1))/x(i)^c;
        z(j) = z(j) + buf * f_i/f_j;
        A(m, 2) = trapz(x, z./x.^c);
    end
    assert(all(abs(A(:, 1) - A(:, 2)) < 1e-8), 'Check reassignment factor');

    figure(1);
    clf; cla;
    plot(x, y);
    hold on
    plot(x, z, '-.m');
    stem(x(1), y(1), 'Color', "#0072BD", 'Marker', 'none');
    stem(x(end), y(end), 'Color', "#0072BD", 'Marker', 'none');
    hold off
    title(sprintf('Area = %3.4f, %3.4f', A(M, 1), A(M, 2)));
end