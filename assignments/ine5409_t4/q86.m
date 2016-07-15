% error = 1e-6
% f(x) = erf(x)

const = 2 / sqrt(pi);
x = 1; a = 0; b = 1;
exact = erf(x);

printf('8.6a) trapezoidal rule')

% error = max(abs(f^(2)(x))_{[a, b]}) * h^2 * (b - a) * 1/12

% f^(2)(x) = (-4 * x * e^(-x^2)) / sqrt(pi)
% max(abs(f^(2)(0)), abs(f^(2)(1))) ≈ 0.8302149948

% 1e-6 ≈ 0.8302149948 * h^2 * (1 - 0) * 1/12
% h ≈ sqrt(1e-6 * 12 / 0.8302149948)
% h ≈ 0.0038018531

% h = (b - a) / n
% n ≈ (1 - 0) / 0.0038018531
% n ≈ 263.0296236

nT = 264;
Tn = trapezoidal_rule(nT, a, b) * const;
Tnn = trapezoidal_rule(2 * nT, a, b) * const;
T_exact_error = abs(Tn - exact);
T_est_error = abs(Tn - Tnn);

printf('\n  number of points n      = %d', nT)
printf('\n  exact value of erf(1)   = %.20f', exact)
printf('\n  value with n points     = %.20f', Tn)
printf('\n  value with 2n points    = %.20f', Tnn)
printf('\n  maximum exact error     = %.20f', T_exact_error)
printf('\n  maximum estimated error = %.20f', T_est_error)

printf('\n\n8.6b) Simpson`s rule')

% error = max(abs(f^(4)(x)_{[a, b]})) * h^4 * (b - a) * 1/180

% f^(4)(x) = (-8 * x * e^(-x^2) * (2x^2 - 3)) / sqrt(pi)
% max(abs(f^(4)(0)), abs(f^(4)(1))) ≈ 1.6604299896

% 1e-6 ≈ 1.6604299896 * h^4 * (1 - 0) * 1/180
% h ≈ (1e-6 * 180 / 1.6604299896)^(1/4)
% h ≈ 0.1020382458

% h = (b - a) / n
% n ≈ (1 - 0) / 0.1020382458
% n ≈ 9.8002468726

nS = 10;
Sn = simpson_rule(nS, a, b) * const;
Snn = simpson_rule(2 * nS, a, b) * const;
S_exact_error = abs(Sn - exact);
S_est_error = abs(Sn - Snn);

printf('\n  number of points n      = %d', nS)
printf('\n  exact value of erf(1)   = %.20f', exact)
printf('\n  value with n points     = %.20f', Sn)
printf('\n  value with 2n points    = %.20f', Snn)
printf('\n  maximum exact error     = %.20f', S_exact_error)
printf('\n  maximum estimated error = %.20f', S_est_error)

printf('\n\n8.6c) Gauss-Legendre quadrature')

% error ≤ max(abs(f^(2m)(x))_{[a, b]}) * (b - a)^(2*m - 1)
%           * (m!)^4 / ((2*m + 1) * [(2*m)!]^3)
% (special case of ω(x) = 1)

% m = 1
% f^(2)(x) = (-4 * x * e^(-x^2)) / sqrt(pi)
% max(abs(f^(2)(0)), abs(f^(2)(1))) ≈ 0.8302149948
% error ≈ 0.8302149948 * (1 - 0)^(2*1 - 1) * (1!)^4 / ((2*1 + 1) * [(2*1)!]^3)
%       ≈ 0.8302149948 * 1 * 1 / (3 * 2^3)
%       ≈ 0.8302149948 / 24
%       ≈ 0.0345922914

% m = 2
% f^(4)(x) = (-8 * x * e^(-x^2) * (2x^2 - 3)) / sqrt(pi)
% max(abs(f^(4)(0)), abs(f^(4)(1))) ≈ 1.6604299896
% error ≈ 1.6604299896 * (1 - 0)^(2*2 - 1) * (2!)^4 / ((2*2 + 1) * [(2*2)!]^3)
%       ≈ 1.6604299896 * 1 * 16 / (5 * 24^3)
%       ≈ 1.6604299896 * 16 / 69120
%       ≈ 0.0003843587

% m = 3
% f^(6)(x) = (-16 * x * e^(-x^2) * (4*x^4 - 20*x^2 + 15)) / sqrt(pi)
% max(abs(f^(6)(0)), abs(f^(6)(1))) ≈ 3.3208599793
% error ≈ 3.3208599793 * (1 - 0)^(2*3 - 1) * (3!)^4 / ((2*3 + 1) * [(2*3)!]^3)
%       ≈ 3.3208599793 * 1 * 1296 / (7 * 720^3)
%       ≈ 3.3208599793 * 1296 / 2612736000
%       ≈ 0.0000016472

% m = 4
% f^(8)(x) = (-32 * x * e^(-x^2) * (8*x^6 - 84*x^4 + 210*x^2 - 105)) / sqrt(pi)
% max(abs(f^(8)(0)), abs(f^(8)(1))) ≈ 192.6098788
% error ≈ 192.6098788 * (1 - 0)^(2*4 - 1) * (4!)^4 / ((2*4 + 1) * [(2*4)!]^3)
%       ≈ 192.6098788 * 1 * 331776 / (9 * 40320^3)
%       ≈ 192.6098788 * 331776 / 589934886912000
%       ≈ 0.0000001083

m = 4;
[Gm, c, xm, ym] = gauss_legendre(m, a, b);
Gm *= const;
Gmm = gauss_legendre(m + 1, a, b) * const;
G_exact_error = abs(Gm - exact);
G_est_error = abs(Gm - Gmm);

printf('\n  number of points m      = %d', m)
printf('\n  exact value of erf(1)   = %.20f', exact)
printf('\n  value with m points     = %.20f', Gm)
printf('\n  value with m + 1 points = %.20f', Gmm)
printf('\n  maximum exact error     = %.20f', G_exact_error)
printf('\n  maximum estimated error = %.20f', G_est_error)

printf('\n\n8.6d) interpolating polynomial from Gauss-Legendre where n = m - 1')

printf('\n  P_n(x) = ')
for i = 1 : m - 1
    printf('%.9f*x^%d + ' , c(i), m - i)
end
printf('%.9f', c(m))

printf('\n\n8.6e) integral of the interpolating polynomial P_n(x)')

% polyint takes an array with the coefficients of the polynomial
integral = polyint(c);
% polyval applies a value to a polynomial
area_value = polyval(integral, b) - polyval(integral, a);
% should be equal to Gm, and indeed, it is
printf('\n  area of interval [a, b] = %.20f', area_value)

printf('\n\n8.6f) plot of erf(x) and P_n(x)')

x = a : 0.01 : b;
y = const .* exp(-x .^ 2);
yp = polyval(c, x);

%plot(
%    x,  y,  '-k;erf(x);', 'linewidth', 4,
%    x,  yp, '-r;P_n(x);', 'linewidth', 2,
%    xm, ym, 'bx;Gauss-Legendre points;', 'linewidth', 2
%);

printf('\n\n')
