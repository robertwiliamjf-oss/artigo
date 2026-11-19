% =========================================================================
% ANALISE DE CARREGAMENTO COM CONTROLE MULTIDIMENSIONAL DO SVC
% Varredura temporal (24 horas) com fatores de desequilíbrio variáveis.
% Versão estendida: simula múltiplos casos e salva dados completos.
% =========================================================================
clc; clear; close all;
t_ini = tic;
lambda_val=1;
global lii SVC THDv V alpha lambda
global FAKTOR Ton1s Ton2s
global des1 des2 des3

% ========== PARÂMETROS DO SISTEMA ==========
n_barras = 33;
alpha_inicial = 110;
fatores_iniciais = [0.5,0.5,0.5];
lambda_vector = 1;          % apenas um valor de lambda (nominal)
n_pontos = length(lambda_vector);

% ========== DEFINIÇÃO DOS FATORES DE DESEQUILÍBRIO DIÁRIOS ==========
des1_vec = [1.05, 1.08, 1.10, 1.12, 1.10, 1.08, ...
            1.05, 1.02, 0.98, 0.95, 0.93, 0.92, ...
            0.93, 0.95, 0.98, 1.00, 1.02, 1.03, ...
            1.02, 1.00, 0.98, 0.96, 0.95, 0.94];
des2_vec = [0.95, 0.94, 0.93, 0.92, 0.93, 0.95, ...
            0.98, 1.02, 1.05, 1.08, 1.10, 1.12, ...
            1.10, 1.08, 1.05, 1.02, 0.98, 0.95, ...
            0.93, 0.92, 0.93, 0.95, 0.97, 0.98];
des3_vec = [0.92, 0.91, 0.90, 0.89, 0.88, 0.90, ...
            0.93, 0.96, 0.98, 1.00, 1.02, 1.04, ...
            1.06, 1.08, 1.10, 1.12, 1.15, 1.18, ...
            1.16, 1.12, 1.08, 1.02, 0.98, 0.94];
combinacoes = [des1_vec(:), des2_vec(:), des3_vec(:)];
n_horas = size(combinacoes, 1);
fprintf('Total de horas simuladas: %d\n', n_horas);

% ========== CONFIGURAÇÃO DO OTIMIZADOR ==========
options = optimoptions('fmincon', 'Display', 'none', ...
                       'Algorithm', 'sqp', ...
                       'MaxIterations', 400, ...
                       'TolFun', 1e-8, ...
                       'TolX', 1e-8, ...
                       'MaxFunctionEvaluations', 1e5);

% ========== CASOS A SEREM SIMULADOS ==========
% Definição dos casos
casos = struct();
% Caso 0: sem SVC (apenas cargas, sem controle)
casos(1).nome = 'sem_svc';
casos(1).SVC = 0;
casos(1).otimizar = false;
% Caso 1: apenas capacitores fixos (mesmo que SVC=0, mas mantemos separado)
casos(2).nome = 'fixo_capacitor';
casos(2).SVC = 0;
casos(2).otimizar = false;
% Caso 2: ângulo fixo (100°)
casos(3).nome = 'angulo_fixo';
casos(3).SVC = 1;
casos(3).otimizar = false;
casos(3).alpha_fixo = 100;
casos(3).fatores_fixo = [1,1,1];
% Caso 3: controle balanceado (Ângulos iguais nas três fases, otimizados)
casos(4).nome = 'balanceado';
casos(4).SVC = 1;
casos(4).otimizar = true;
casos(4).usar_fatores = false;   % fixa FAKTOR = [1,1,1]
% Caso 4: proposto (com fatores de fase otimizados)
casos(5).nome = 'proposto';
casos(5).SVC = 1;
casos(5).otimizar = true;
casos(5).usar_fatores = true;

n_casos = length(casos);

% ========== PRE‑ALOCAÇÃO DE DADOS (estrutura completa) ==========
dados_finais = struct();
dados_finais.n_barras = n_barras;
dados_finais.n_horas = n_horas;
dados_finais.lambda_vec = lambda_vector;
dados_finais.combinacoes = combinacoes;

% Inicializa campos para cada caso
for c = 1:n_casos
    nome = casos(c).nome;
    dados_finais.(nome).tensoes   = zeros(n_barras, 3, n_horas);
    dados_finais.(nome).thdv      = zeros(n_barras, 3, n_horas);
    dados_finais.(nome).vuf       = zeros(2, n_horas);   % barras 17 e 32
    dados_finais.(nome).alpha1    = zeros(1, n_horas);
    dados_finais.(nome).alpha2    = zeros(1, n_horas);
    dados_finais.(nome).fatores   = zeros(3, n_horas);
end

% ========== INICIALIZAÇÃO COMUM (caso base neutro para V_ini) ==========
des1 = 1; des2 = 1; des3 = 1;
lambda=1;
alpha1=alpha_inicial;
alpha2=alpha1;
Header_1;
V_ini = V;

% ========== LOOP PRINCIPAL SOBRE CASOS ==========
for c = 1:n_casos
    caso = casos(c);
    fprintf('\n--- Simulando caso: %s ---\n', caso.nome);
    
    % Pre‑alocação local
    tensoes_temp = zeros(n_barras, 3, n_horas);
    thdv_temp    = zeros(n_barras, 3, n_horas);
    vuf_temp     = zeros(2, n_horas);
    alpha1_temp  = zeros(1, n_horas);
    alpha2_temp  = zeros(1, n_horas);
    fatores_temp = zeros(3, n_horas);
    
    % Configuração específica do caso
    SVC = caso.SVC;
    
    for t = 1:n_horas
        des1 = combinacoes(t,1);
        des2 = combinacoes(t,2);
        des3 = combinacoes(t,3);
        
        if ~caso.otimizar
            % Sem otimização – usa valores fixos
            if strcmp(caso.nome, 'angulo_fixo')
                alpha1 = caso.alpha_fixo;
                alpha2 = caso.alpha_fixo;
                FAKTOR = caso.fatores_fixo;
                Ton1s = alpha1;
                Ton2s = alpha2;
            else
                % Sem SVC ou apenas capacitores: ângulos não influenciam
                alpha1 = 90;
                alpha2 = 90;
                FAKTOR = [1,1,1];
                Ton1s = alpha1;
                Ton2s = alpha2;
            end
            Header_1;
            
            % Armazenamento
            alpha1_temp(t) = alpha1;
            alpha2_temp(t) = alpha2;
            fatores_temp(:,t) = FAKTOR;
            for barra = 1:n_barras
                tensoes_temp(barra, :, t) = [V(barra,3), V(barra,4), V(barra,5)];
                thdv_temp(barra, :, t)   = [THDv(barra,1), THDv(barra,2), THDv(barra,3)];
            end
            [VUF, ~] = calcular_vuf(V, T, [17,32]);
            vuf_temp(:,t) = VUF;
        else
            % Com otimização
            % Reinicia ângulos a cada hora

            FAKTOR = fatores_iniciais;
            Ton1s = alpha1;
            Ton2s = alpha2;
            
            % Define se usa fatores
            if ~caso.usar_fatores
                % Controle balanceado: fixa FAKTOR = [1,1,1] e otimiza apenas alphas
                % Mas o fmincon precisa receber 7 variáveis; podemos fixar os fatores
                % e otimizar apenas as 4 primeiras.
                % Abordagem mais simples: executar fmincon com 7 variáveis, mas com
                % fatores fixos (limites lb(5:7)=ub(5:7)=1) e remover da otimização.
                lb = [1, 0, 90, 90, 1, 1, 1];
                ub = [1, 0, 180, 180, 1, 1, 1];
                % Ponto inicial inclui fatores=1
                x0 = [1, 0, alpha1, alpha2, 1, 1, 1];
                FAKTOR_prev = [1,1,1];
                % Função objetivo sem argumento extra para fatores? A fob já usa FAKTOR_prev.
                [x, fval, exitflag] = fmincon(@(X) fob(X, SVC, V, V_ini, alpha1, alpha2, des1, des2, des3, FAKTOR_prev), ...
                                              x0, [], [], [], [], lb, ub, ...
                                              @(X) in(X, V, THDv, lii, SVC, Ta, Tb, Tc, des1, des2, des3), options);
            else
                % Controle proposto: otimiza todos os 7 parâmetros
                x0 = [lambda_val, 0, alpha1, alpha2, FAKTOR];
                lb = [lambda_val, 0, 90, 90, 0, 0, 0];
                ub = [lambda_val, 0, 180, 180, 100, 100, 100];
                FAKTOR_prev = [1,1,1];
                [x, fval, exitflag] = fmincon(@(X) fob(X, SVC, V, V_ini, alpha1, alpha2, des1, des2, des3, FAKTOR_prev), ...
                                              x0, [], [], [], [], lb, ub, ...
                                              @(X) in(X, V, THDv, lii, SVC, Ta, Tb, Tc, des1, des2, des3), options);
            end
            
            if exitflag >= 0
                lambda = x(1);
                lbd = x(2);
                alpha1 = x(3);
                alpha2 = x(4);
                if caso.usar_fatores
                    FAKTOR = x(5:7);
                else
                    FAKTOR = [1,1,1];
                end
                Ton1s = alpha1;
                Ton2s = alpha2;
            else
                fprintf('Falha na hora %02d:00 (exitflag=%d).\n', t-1, exitflag);
            end
            
            Header_1;
            
            % Armazenamento
            alpha1_temp(t) = alpha1;
            alpha2_temp(t) = alpha2;
            fatores_temp(:,t) = FAKTOR;
            for barra = 1:n_barras
                tensoes_temp(barra, :, t) = [V(barra,3), V(barra,4), V(barra,5)];
                thdv_temp(barra, :, t)   = [THDv(barra,1), THDv(barra,2), THDv(barra,3)];
            end
            [VUF, ~] = calcular_vuf(V, T, [17,32]);
            vuf_temp(:,t) = VUF;
        end
    end
    
    % Salva dados do caso na estrutura final
    dados_finais.(caso.nome).tensoes = tensoes_temp;
    dados_finais.(caso.nome).thdv    = thdv_temp;
    dados_finais.(caso.nome).vuf     = vuf_temp;
    dados_finais.(caso.nome).alpha1  = alpha1_temp;
    dados_finais.(caso.nome).alpha2  = alpha2_temp;
    dados_finais.(caso.nome).fatores = fatores_temp;
end

% ========== SALVAR DADOS COMPLETOS ==========
save('dados_completos_artigo.mat', '-struct', 'dados_finais');
fprintf('Dados salvos em dados_completos_artigo.mat\n');
fprintf('Tempo total: %.2f minutos\n', toc(t_ini)/60);