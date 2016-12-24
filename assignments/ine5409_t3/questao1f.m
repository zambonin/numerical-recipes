function x = questao1f()

	x = -1 : 0.01 : +1;
	exato = sin(x);

	approxCheb = (x .* (383/384)) - ((x .^ 3) .* (5/32));
    erroCheb = abs(approxCheb .- exato);

    approxPade = (x .* (5121/5120)) - ((x .^ 3) .* (959/5760)) + ((x .^ 5) .* (23/2880));
    erroPade = abs(approxPade .- exato);

	plot(exato, erroCheb, "b; Chebyshev;",
		 exato, erroPade, "g; Pade;")

end
