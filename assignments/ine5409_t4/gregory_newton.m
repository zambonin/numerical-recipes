function yPn = gregory_newton(p, n, x, y)
    
    k = 1;
    for i = 1 : n
        DDy(i, k) = (y(i+1) - y(i)) / (x(i+1) - x(i));
    end
    for k = 2 : n
        for i = 1 : (n-k+1)
            DDy(i, k) = (DDy(i+1, k-1) - DDy(i, k-1)) / (x(i+k) - x(i));
        end
    end

    for i = 1 : length(p)
        t_sum = y(1);
        for k = 1 : n
            t_sum += DDy(1, k) * prod(p(i) - x(1:k));
        end
        yPn(i) = t_sum;
    end

end
