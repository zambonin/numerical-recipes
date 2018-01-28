function [a b c d S]=fSplinea(m,x,y,h)
    % teremos 'm' splines cúbicas, uma para cada intervalo, formando um sistema
    % tridiagonal de 'm-1' equacoes para S (derivadas de segunda ordem)
    % Splines naturais/lineares NAS PONTAS da figura S1=0 e Sm+1=0
    %1a. eq. i=2
    t(2)=0;
    r(2)=2*h(1)+2*h(2);
    d(2)=h(2);
    b(2)=6*((y(3)-y(2))/h(2)-(y(2)-y(1))/h(1)); %para S1
    for i=3:m-1
    	t(i)=h(i-1);
    	r(i)=2*(h(i-1)+h(i));
    	d(i)=h(i);
    	b(i)=6*((y(i+1)-y(i))/h(i)-(y(i)-y(i-1))/h(i-1)); %para Si
    end
    %ultima eq. i=m
    t(m)=h(m-1);
    r(m)=(2*h(m-1)+2*h(m));
    d(m)=0;
    b(m)=6*((y(m+1)-y(m))/h(m)-(y(m)-y(m-1))/h(m-1)); %para Sm
    %Sistemas de eqs. tridiagonal com eqs. de i=2 a m
    for i=3:m
        aux=t(i)/r(i-1);t(i)=0;
        r(i)=r(i)-aux*d(i-1);
        b(i)=b(i)-aux*b(i-1);
    end
    S(m)=b(m)/r(m);
    for i=m-1:-1:2
      S(i)=(b(i)-d(i)*S(i+1))/r(i);
    end
    % Splines naturais/lineares NAS PONTAS da figura S1=0 e Sm+1=0
    S(1)=0; S(m+1)=0;
    S;
    % Calcula-se cada conjunto de coeficientes a, b, c, d dos polinimois de 3o. grau
    for i=1:m
         a(i)=(S(i+1)-S(i))/(6*h(i));
         b(i)= S(i)/2;
         c(i)=(y(i+1)-y(i))/h(i)-(S(i+1)+2*S(i))*h(i)/6;
         d(i)= y(i);
    end
end