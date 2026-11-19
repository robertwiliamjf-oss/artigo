% =========================================================================
% ANALISE DE CARREGAMENTO COM CONTROLE MULTIDIMENSIONAL DO SVC
% Objetivo: minimizar angulos de disparo, mantendo tensoes e THD nos limites.
% =========================================================================
clc; clear; close all;
t_ini = tic;

global lii SVC THDv V alpha lambda

% Parametros do sistema
n_barras = 33;          % numero de barras
alpha_inicial = 90;     % angulo inicial (ponto de repouso)
alpha1_prev = 90;
alpha2_prev = 90;
% % ========== 1. DEFINICAO DO VETOR DE CARGA (curva simples) ==========
% Rampa de subida e descida (evita complexidade inicial)
% Par�metros
passo = 0.05;
lambda_min = 1;
lambda_max = 2;
pontos_plato = 5;   % quantas vezes repetir o valor 2

% Vetores
subida = lambda_min:passo:lambda_max;          % 1, 1.1, ..., 2
plato = lambda_max * ones(1, pontos_plato);    % [2,2,2,2,2]
descida = (lambda_max-passo):-passo:lambda_min; % 1.9, 1.8, ..., 1

% % Vetor final
lambda_vector = [subida,plato, descida];
lambda_vector = curva_de_carga(1.0, 1.65);
% lambda_vector = gerar_lambda_duck(100, 0.9, 2, ...
%     'trough_frac', 0.3, ...        % vale baixo (0.9 + 0.3*(1.1) = 1.23)
%     'p1', 0.2, ...                 % primeiro pico em 20% do tempo
%     'p2', 0.5, ...                 % vale no meio do per�odo (mais distante do pico)
%     'p3', 0.8, ...                 % segundo pico em 80%
%     'smooth', 'spline', ...
%     'waviness', 0.05, ...          % ondula��o pequena (5% da amplitude)
%     'waves', 2);
n_pontos = length(lambda_vector);

% ========== 2. PRE-ALOCACAO DE DADOS ==========
tensoes_3d = zeros(n_barras, 3, n_pontos);   % (barra, fase, ponto)
thdv_3d   = zeros(n_barras, 3, n_pontos);    % (barra, fase, ponto)
lambdas_armazenados = zeros(1, n_pontos);
LLLA = zeros(1, n_pontos);
ALPHA_1 = zeros(1, n_pontos);
ALPHA_2 = zeros(1, n_pontos);
lbd_values = zeros(1, n_pontos);

% ========== 3. INICIALIZA��ES ==========
lambda = 1;
alpha1 = alpha_inicial;
alpha2 = alpha_inicial;
SVC = 1;
Header_1;
V_ini = V;

lambda_anterior = lambda_vector(1);   % para detectar descida (opcional)

% ========== 4. CONFIGURA��O DO OTIMIZADOR ==========
options = optimoptions('fmincon', 'Display', 'iter', ...
                       'Algorithm', 'sqp', ...
                       'MaxIterations', 20, ...
                       'TolFun', 1e-6, ...
                       'TolX', 1e-6);

% ========== 5. LOOP PRINCIPAL ==========
h_waitbar = waitbar(0, 'Processando pontos de carga...');
for count = 1:n_pontos
    lambda_val = lambda_vector(count);
    waitbar(count/n_pontos, h_waitbar, ...
            sprintf('Ponto %d de %d (lambda = %.3f)', count, n_pontos, lambda_val));
    
    % (Opcional) Se a carga est� diminuindo e os �ngulos est�o muito altos,
    % podemos reiniciar para facilitar a descida. Mas com o novo objetivo,
    % isso pode n�o ser necess�rio.
    % if lambda_val < lambda_anterior && (alpha1 > 150 || alpha2 > 150)
    %     alpha1 = alpha_inicial;
    %     alpha2 = alpha_inicial;
    %     % fprintf('Reinicializando angulos em lambda = %.3f\n', lambda_val);
    % end
    % 
    % Chama a otimiza��o
    x0 = [lambda_val 0 alpha1 alpha2];
   
    lb = [lambda_val 0 90 90];
    ub = [lambda_val 0 90 90];
 [x, fval, exitflag] = fmincon(@(X) fob(X, SVC, V, V_ini, alpha1_prev, alpha2_prev), ...
                                  x0, [], [], [], [], lb, ub, ...
                                  @(X) in(X, V, THDv, lii, SVC), options);
    if exitflag <= 0
    fprintf('Aviso: falha em lambda=%.3f (exitflag=%d). Mantendo �ngulos anteriores.\n', lambda_val, exitflag);
    % N�o atualiza alpha1, alpha2 � eles continuam como estavam
else
    alpha1 = x(3);
    alpha2 = x(4);
end
    % Atualiza vari�veis
    lambda = x(1);
    lbd = x(2);
    alpha1_prev = alpha1;
    alpha2_prev = alpha2;
    Header_1;   % Atualiza as vari�veis globais V e THDv com os �ngulos �timos
    
    % Armazena resultados
    LLLA(count) = lambda;
    lbd_values(count) = lbd;
    ALPHA_1(count) = alpha1;
    ALPHA_2(count) = alpha2;
    lambdas_armazenados(count) = lambda;
    
    % Salva tensoes e THDv das 3 fases para todas as barras
    for barra = 1:n_barras
        tensoes_3d(barra, :, count) = [V(barra, 3), V(barra, 4), V(barra, 5)];
        thdv_3d(barra, :, count)   = [THDv(barra, 1), THDv(barra, 2), THDv(barra, 3)];
    end
    
    lambda_anterior = lambda_val;
end
close(h_waitbar);

% ========== 6. SALVAR DADOS ==========
save('dados_analise.mat', 'LLLA', 'ALPHA_1', 'ALPHA_2', 'lbd_values', ...
     'tensoes_3d', 'thdv_3d', 'lambdas_armazenados', 'n_barras');

fprintf('Dados salvos em dados_analise.mat\n');
fprintf('Tempo total: %.2f minutos\n', toc(t_ini)/60);

% ========== 7. PLOTAR RESULTADOS (chamar fun��o separada) ==========
plotar_graficos_selecionados();   % fun��o definida anteriormente