printf('7.3a) functions for temperature x volume')

m = 7;
T = [13.9 37.0 67.8 79.0 85.5 93.1 99.2];
V = [1.04 1.18 1.29 1.35 1.28 1.21 1.06];

xp = T(1) : (T(m) - T(1)) / 100 : T(m);

a1 = poly_adjust(1, m, T, V);
a2 = poly_adjust(2, m, T, V);

printf('\n  P_1(x) = %.15f*x^1 + %.15f', a1(1), a1(2))
printf('\n  P_2(x) = %.15f*x^2 + %.15f*x^1 + %.15f', a2(1), a2(2), a2(3))

printf('\n\n7.3b) interpolating polynomial')

Pn = gregory_newton(xp, m - 1, T, V);

printf('\n  P_n(x) = ')
for i = 1 : floor(m/2)
    printf('%.5f*x^%d + ' , Pn(i), m - i + 1)
end
printf('\n           ')
for i = ceil(m/2) : m
    printf('%.5f*x^%d + ' , Pn(i), m - i + 1)
end
printf('%.5f', Pn(m+1))

printf('\n\n7.3c) plot with points and functions')

y1 = horner(1, poly_adjust(1, m, T, V), xp);
y2 = horner(2, poly_adjust(2, m, T, V), xp);

%plot(
%    xp, y1, '-k;P_1(T);',    'linewidth', 2,
%    xp, y2, '-r;P_2(T);',    'linewidth', 2,
%    xp, Pn, '-g;P_6(T);',    'linewidth', 2,
%    T,  V,  'bx;real data;', 'linewidth', 2
%)

printf('\n\n7.3d) The linear adjustment does not account for the small')
printf('\n  curvature between the points, and the polynomial with the largest')
printf('\n  order has no pattern at all. The best choice is the parabolic')
printf('\n  adjustment.')

printf('\n\n')
