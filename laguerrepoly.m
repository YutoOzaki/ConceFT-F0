function y = laguerrepoly(x, k, c)
    if k >= 0
        m = (0:k)';
        C = (-1).^m.*gamma(k + c + 1)./(gamma(c + m + 1).*gamma(k - m + 1))./factorial(m);
        y = sum(C .* x.^m, 1);
    else
        y = 0.*x;
    end
end