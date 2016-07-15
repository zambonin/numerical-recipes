printf('7.1a) V(T) = ln(a + b * T^2)')

m = 6;
T = [0.20 0.40 0.60 0.80 0.90 1.00];
V = [0.04 0.14 0.30 0.45 0.61 0.69];

xi = [0.01 0.01];

s = newton_raph_nl(xi);
printf('\n  a = %.20f', s(1))
printf('\n  b = %.20f', s(2))

printf('\n\n7.1b) standard deviations')

d = abs(log(s(1) + s(2) .* T.^2) .- V);

printf('\n  d =')
for i = 1 : length(d)
    printf(' %.9f', d(i))
end

printf('\n\n7.1c) plot with points and functions')

xp = T(1) : (T(m) - T(1)) / 100 : T(m);

f = log(s(1) + s(2) .* xp.^2);
Pn = gregory_newton(xp, m - 1, T, V);

% looking at the plot, it should be clear that the adjusted function was better
%plot(
%    xp, f,  '-k;V(T);',      'linewidth', 2,
%    xp, Pn, '-r;P_n(T);',    'linewidth', 2,
%    T,  V,  'bx;real data;', 'linewidth', 2
%)

printf('\n\n')
