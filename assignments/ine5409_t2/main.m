clc
format long

poly = transpose([
    complex(1, 0),
    complex(-5, -1),
    complex(9.99, 5),
    complex(-9.97, -9.99),
    complex(4.97, 9.97),
    complex(-0.99, -4.97),
    complex(0., 0.99)
]);

[x M it] = orig_roots(poly);

printf("\n    b)  Os valores para o limite mínimo dos restos testados foram\n")
printf("        os seguintes:\n")
printf("        r_limite\titerações\tnúmero de raízes\tmultiplicidade\n")
printf("          1e-1  \t   68    \t       2        \t  4, 4       \n")
printf("          1e-2  \t   72    \t       2        \t  4, 3       \n")
printf("          1e-3  \t   138   \t       3        \t  3, 2, 1    \n")
printf("          1e-4  \t   207   \t       4        \t  3, 1, 1, 1 \n")
printf("          1e-5  \t   218   \t       4        \t  3, 1, 1, 1 \n")
printf("          1e-6  \t   244   \t       4        \t  3, 1, 1, 1 \n")

printf("\n    f)  Número total de iterações efetuadas: "), it

printf("\n    g)  Raízes refinadas convergidas e multiplicidades associadas: ")
own = transpose(x), M

printf("    h)  Raízes do Octave e Wolfram Alpha: ")
octave = roots(poly)
wolfram = [
    0.9,
    0.999959,
    1.1,
    complex(0, 1),
    complex(1.00002, -0.0000351069),
    complex(1.00002, +0.0000351069)
]

printf("    i)  Descartando a precisão absurda do Wolfram Alpha, é possível\n")
printf("        verificar que as raízes calculadas por este programa são   \n")
printf("        mais precisas -- é aparente que o número de casas decimais \n")
printf("        válidas após à vírgula é maior no primeiro resultado. Isso \n")
printf("        acontece em virtude das diversas estratégias de aproximação\n")
printf("        por diferentes cotas, refinamentos e redução de grau do    \n")
printf("        polinômio original.\n")
