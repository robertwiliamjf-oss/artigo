function J = fob(X, SVC, V, V_ini, alpha1_prev, alpha2_prev)
    lambda = X(1);
    lbd = X(2);
    alpha1 = X(3);
    alpha2 = X(4);
    
    Header_1;
    
    desvio = 0;
    for b = [17, 32]
        for f = 3:5
            V_target = V_ini(b, f);   % valor inicial para esta barra/fase
            desvio = desvio + (V(b, f) - V_target)^2;
        end
    end
    
    custo_angulo = (alpha1^2 + alpha2^2) / 1e3;
    J = 1e6 * desvio + custo_angulo;
end