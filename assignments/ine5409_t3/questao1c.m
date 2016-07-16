function [C erro] = questao1c()

    printf("\n1c)\n");

    % O resto da série de Maclaurin, dado por R_n(x), pode ser calculado com a
    % fórmula abaixo, e deve respeitar a ordem de O(10^-3). Então,

    %     R_n(x) = (max(abs(f^(n+1)(x))*(x - b)^(n+1))) / (n + 1)! (x em [-1, 1])
    %     R_n(x) <= sqrt(10) * 10^-2

    % Aplicando n = 5 por tentativas, tem-se
        
    %     R_n(x)  = (max(abs(sen(x))*(1 - 0)^6)) / 6!
    %             = sen(1) / 6!
    %             = 0.0011687097...

    % Tal valor respeita a inequação acima, então truncar a série de Maclaurin
    % no quinto termo não afetará o resultado para esta tolerância.

    % Então, seja a série de Maclaurin truncada no quinto termo para sen(x)

    %     M_x = x - (x^3 / 3!) + (x^5 / 5!)

    % e as potências de x escritas em função de polinômios de Chebyshev (e vice-versa)

    %     x = T_1     x^3 = (T_3 + 3T_1) / 4      x^5 = (T_5 + 5T_3 + 10T_1) / 16

    %     T_1 = 2x^2 - 1      T_3 = 4x^3 - 3x     T_5 = 16x^5 - 20x^3 + 5x

    % Substituindo estas potências, tem-se

    %     C_T  = T_1 - (((T_3 + 3T_1) / 4) / 3!) + (((T_5 + 5T_3 + 10T_1) / 16) / 5!)
    %          = T_1 - ((T_3 + 3T_1) / 24) + ((T_5 + 5T_3 + 10T_1) / 1920)
    %          = T_1 - (T_3/24) - (T_1/8) + (T_5/1920) + (T_3/384) + (T_1/192)
    %          = T_1 (1 - 1/8 + 1/192) + T_3 (1/384 - 1/24) + T_5 (1/1920)
    %          = T_1 (169/192) + T_3 (-5/128) + T_5 (1/1920)
        
    % Como 1/1920 < O(10^-3), então T_5 pode ser truncado.

    %     C_T = T_1 (169/192) + T_3 (-5/128)

    % Substituindo os polinômios de Chebyshev para a variável original x, tem-se

    %     C_x = x (169/192) + (4x^3 - 3x) (-5/128)
    %         = x (169/192) - x^3 (20/128) + x (15/128)
    %         = x (383/384) - x^3 (20/128)
    %         = x (383/384) - x^3 (5/32)

    x = [-1, -0.5, 0, 0.5, 1]
    exato = sin(x);
    approx = (x .* (383/384)) - ((x .^ 3) .* (5/32));

    C = [383/384, 5/32];
    erro = max(abs(approx .- exato));

end
