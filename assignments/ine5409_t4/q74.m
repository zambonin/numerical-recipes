printf('7.4a) V(T) = a * T + b * cos(T)')

m = 4;
T = [0.00 0.39 0.78 1.18];
V = [0.99 0.92 0.71 0.38];

A = zeros(2, 3);

for i = 1 : m
	A(1,1) += T(i)^2;
	A(1,2) += (T(i) * cos(T(i)));
	A(1,3) += (T(i) * V(i));
	A(2,2) += cos(T(i))^2;
	A(2,3) += (V(i) * cos(T(i)));
end
A(2,1) = A(1,2);

a = gauss_elimination(2, A);
printf('\n  a = %.20f', a(1))
printf('\n  b = %.20f', a(2))
 
printf('\n\n7.4b) histogram of local standard deviations')

d = abs(a(1) .* T .+ a(2) .* cos(T) .- V);

printf('\n  d =')
for i = 1 : length(d)
    printf(' %.10f', d(i))
end

%bar(T, d);

printf('\n\n')
