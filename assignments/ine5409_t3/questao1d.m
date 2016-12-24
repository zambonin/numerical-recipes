function [C erro] = questao1d()

    printf("\n1d)\n");

    % Seja a série de Maclaurin truncada no quinto termo para sen(x)

    %     M_x = x - (x^3 / 3!) + (x^5 / 5!) - (x^7 / 7!)

    % e as potências de x escritas em função de polinômios de Chebyshev (e vice-versa)

    %   x = T_1
    %   x^3 = (T_3 + 3T_1) / 4
    %   x^5 = (T_5 + 5T_3 + 10T_1) / 16
    %   x^7 = (T_7 + 7T_5 + 21T_3 + 35T_1) / 64

    %   T_1 = x
    %   T_3 = 4x^3 - 3x
    %   T_5 = 16x^5 - 20x^3 + 5x
    %   T_7 = 64x^7 - 112x^5 + 56x^3 - 7x

    % Substituindo estas potências, tem-se

    %   C_T  = T_1 - (((T_3 + 3T_1) / 4) / 3!) + (((T_5 + 5T_3 + 10T_1) / 16) / 5!)
    %           + (((T_7 + 7T_5 + 21T_3 + 35T_1) / 64) / 7!)
    %        = T_1 - ((T_3 + 3T_1) / 24) + ((T_5 + 5T_3 + 10T_1) / 1920)
    %           - ((T_7 + 7T_5 + 21T_3 + 35T_1) / 322560)
    %        = T_1 (8111/9216) - T_3 (601/15360) + T_5 (23/46080) - T_7 (1/322560)

    % Substituindo os polinômios de Chebyshev para a variável original x, tem-se

    %   C_x = x (8113/9216) - (4x^3 - 3x) (601/15360) + (16x^5 - 20x^3 + 5x) (23/46080)
    %   C_x = x^5 (23/2880) - x^3 (959/5760) + x (5121/5120)

    x = [-1, -0.5, 0, 0.5, 1];
    exato = sin(x);
    approx = (x .* (5121/5120)) - ((x .^ 3) .* (959/5760)) + ((x .^ 5) .* (23/2880));

    C = [5121/5120, 959/5760, 23/2880];
    erro = max(abs(approx .- exato));

end
