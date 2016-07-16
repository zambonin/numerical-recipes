function [a b] = fPade(npade, mpade, c)

    k = npade - mpade;

    for i = 1 : mpade
        for j = 1 : i
            A(i, j) = c(k+i+j);
            A(j, i) = A(i, j);
        end
        A(i, mpade+1) = -c(npade+i+1);
    end

    aux = elimGauss(mpade, A);
    b = fliplr(transpose(aux));
    b(mpade+1 : npade) = 0; 
    b = [1 b];
    a(1) = c(1);

    for i = 2 : npade+1
        S = c(i);
        for j = 1 :i-1
            S = S + b(j+1) * c(i-j);
        end
        a(i) = S;
    end

end
