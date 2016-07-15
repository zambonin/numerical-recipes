function x = newton_raph_nl(xi)

    dx = [1e-6 1e-6];
    tol = 1;
    n_it = 0;

    while (tol > 1e-14 && n_it < 30)
        n_it++;

     	e1 = [xi(1) + dx(1), xi(2)];
        e2 = [xi(1), xi(2) + dx(2)];
        e3 = [xi(1) + dx(1), xi(2)];
        e4 = [xi(1), xi(2) + dx(2)];

		A = [
				(f1a(e1) - f1a(xi)) / dx(1), (f1a(e2) - f1a(xi)) / dx(2), -f1a(xi);
				(f1b(e1) - f1b(xi)) / dx(1), (f1b(e2) - f1b(xi)) / dx(2), -f1b(xi);
			];

		dx = gauss_elimination(2, A);
		x = xi + dx;
		xi = x;
		tol = max(abs(dx));

    end

end
