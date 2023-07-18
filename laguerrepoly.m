function y = laguerrepoly(x, k, c)
    if k >= 0
        m = (0:k)';
        C = (-1).^m.*gamma(k + c + 1)./(gamma(c + m + 1).*gamma(k - m + 1))./factorial(m);
        y = sum(C .* x.^m, 1);
    else
        y = 0.*x;
    end
end

%{
N = 200;
x = rand(1, N);
k = 3;
c = 2;

dy = -laguerrepoly(x, k - 1, c + 1);

eps = 1e-6;
dyhat = (laguerrepoly(x + eps, k, c) - laguerrepoly(x - eps, k, c))./(2*eps);

figure(1);
subplot(2, 1, 1);
plot(dy);
hold on
plot(dyhat, '-.m');
hold off
subplot(2, 1, 2);
histogram(abs(dy - dyhat)./abs(dyhat));
%}