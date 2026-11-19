function [c, ceq] = in(X, V, THDv, lii, SVC, Ta, Tb, Tc,des1,des2,des3)
    % X = [lambda, lbd, alpha1, alpha2]
    lambda  = X(1);
    lbd     = X(2);
    alpha1  = X(3);
    alpha2  = X(4);
    
    % Atualiza o ponto de operação (V, THDv, T e demais variáveis globais)
    Header_1;
    
    % ===== 1. Restrições de tensão =====
    V_min = 0.95;
    V_max = 1.05;
    c_tensao = [ V_min - V(17,3);   V_min - V(17,4);   V_min - V(17,5);
                 V_min - V(32,3);   V_min - V(32,4);   V_min - V(32,5);
                 V(17,3) - V_max;   V(17,4) - V_max;   V(17,5) - V_max;
                 V(32,3) - V_max;   V(32,4) - V_max;   V(32,5) - V_max;
                 90-X(3)*X(5); 90-X(3)*X(6); 90-X(3)*X(7);
                 X(3)*X(5)-180; X(3)*X(6)-180; X(3)*X(7)-180;
                 90-X(4)*X(5); 90-X(4)*X(6); 90-X(4)*X(7);
                 X(4)*X(5)-180; X(4)*X(6)-180; X(4)*X(7)-180;];
    
    % ===== 2. Restrições de THD =====
    thd_limit = 10;   % %
    c_thd = [ THDv(:,1) - thd_limit;
              THDv(:,2) - thd_limit;
              THDv(:,3) - thd_limit ];
    
    % ===== 3. Restrição de desequilíbrio (VUF) =====
    unbalance_limit = 4;   % % (ajuste conforme necessário)
    barras_monitoradas = [17, 32];
    
    % A matriz T (global) contém os ângulos EM RADIANOS (colunas 3,4,5)
    % Extraímos diretamente sem converter
    T_ang = T;                     % ou diretamente T(:,3:5) se preferir
    mags  = V;                     % matriz de magnitudes
    
    % Filtra apenas a frequência fundamental (h=1)
    ind_fund = (mags(:,2) == 1);
    if ~any(ind_fund)
        error('Frequência fundamental (h=1) não encontrada nos dados.');
    end
    
    barra_num = mags(ind_fund, 1);
    Va_mod = mags(ind_fund, 3);
    Vb_mod = mags(ind_fund, 4);
    Vc_mod = mags(ind_fund, 5);
    
    % Ângulos já em radianos: usamos diretamente na exponencial
    Ta_fund = T_ang(ind_fund, 3);
    Tb_fund = T_ang(ind_fund, 4);
    Tc_fund = T_ang(ind_fund, 5);
    
    % Fasores complexos (sem conversão)
    Va_fasorial = Va_mod .* exp(1i * Ta_fund);
    Vb_fasorial = Vb_mod .* exp(1i * Tb_fund);
    Vc_fasorial = Vc_mod .* exp(1i * Tc_fund);
    
    % Operador de Fortescue
    a = exp(1i * 2*pi/3);
    
    % Calcula VUF para cada barra monitorada
    c_unbal = zeros(length(barras_monitoradas), 1);
    for i = 1:length(barras_monitoradas)
        idx = find(barra_num == barras_monitoradas(i), 1);
        if isempty(idx)
            warning('Barra %d não encontrada nos dados da fundamental.', barras_monitoradas(i));
            c_unbal(i) = 0;   % ou trate conforme desejar
            continue;
        end
        
        Va = Va_fasorial(idx);
        Vb = Vb_fasorial(idx);
        Vc = Vc_fasorial(idx);
        
        V_pos = (Va + a*Vb + a^2*Vc) / 3;
        V_neg = (Va + a^2*Vb + a*Vc) / 3;
        
        if abs(V_pos) > 1e-9
            VUF = abs(V_neg) / abs(V_pos) * 100;
        else
            VUF = 0;   % sistema sem tensão positiva → VUF = 0 (ou 100 se preferir)
        end
        
        c_unbal(i) = VUF - unbalance_limit;
    end
    
    % ===== 4. Restrição de resíduo (opcional) =====
    global DI;
    if ~isempty(DI)
        c_res = abs(DI) - 0.001;
    else
        c_res = [];
    end
    
    % ===== Concatena todas as restrições =====
    c = [c_tensao; c_thd; c_unbal];   % pode acrescentar c_res se quiser
    
    ceq = [];
end