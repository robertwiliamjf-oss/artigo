% =========================================================================
% ANALISE DE CARREGAMENTO COM CONTROLE MULTIDIMENSIONAL DO SVC
% Inclui otimiza��o dos fatores de fase (FAKTOR) para �ngulos diferentes por fase.
% =========================================================================
clc; clear; close all;
t_ini = tic;

global lii SVC THDv V alpha lambda
global FAKTOR Ton1s Ton2s   % +++++ adicionado para propagar fatores e �ngulos base

% Parametros do sistema
n_barras = 33;
alpha_inicial = 120;
fatores_iniciais = [0.1,0.1,0.1];   % +++++ fatores iniciais (sem desequil�brio)

% ========== 1. DEFINICAO DO VETOR DE CARGA ==========
passo = 0.05;
lambda_min = 1;
lambda_max = 2;
pontos_plato = 5;
subida = lambda_min:passo:lambda_max;
plato = lambda_max * ones(1, pontos_plato);
descida = (lambda_max-passo):-passo:lambda_min;
lambda_vector = 1.2;
n_pontos = length(lambda_vector);

% ========== 2. PRE-ALOCACAO DE DADOS ==========
tensoes_3d = zeros(n_barras, 3, n_pontos);
thdv_3d   = zeros(n_barras, 3, n_pontos);
lambdas_armazenados = zeros(1, n_pontos);
LLLA = zeros(1, n_pontos);
ALPHA_1 = zeros(1, n_pontos);
ALPHA_2 = zeros(1, n_pontos);
FATORES = zeros(n_pontos, 3);         % +++++ armazenar fatores
lbd_values = zeros(1, n_pontos);

% ========== 3. INICIALIZA��ES ==========
lambda=1;
alpha1 = alpha_inicial;
alpha2 = alpha_inicial;
FAKTOR = fatores_iniciais;
% Ton1s = alpha1;    % +++++ �ngulos base dos SVCs
% Ton2s = alpha2;
SVC = 1;
Header_1;
V_ini = V;

lambda_anterior = lambda_vector(1);

% ========== 4. CONFIGURA��O DO OTIMIZADOR ==========
% +++++ agora com 7 vari�veis
options = optimoptions('fmincon', 'Display', 'iter', ...
                       'Algorithm', 'sqp', ...
                       'MaxIterations', 30, ...
                       'TolFun', 1e-6, ...
                       'TolX', 1e-6);

% ========== 5. LOOP PRINCIPAL ==========
h_waitbar = waitbar(0, 'Processando pontos de carga...');
for count = 1:n_pontos
    lambda_val = lambda_vector(count);
    waitbar(count/n_pontos, h_waitbar, ...
            sprintf('Ponto %d de %d (lambda = %.3f)', count, n_pontos, lambda_val));

    % +++++ Ponto inicial inclui os fatores atuais
    x0 = [lambda_val, 0, alpha1, alpha2, FAKTOR];
    % +++++ Limites: fatores podem variar, ex: 0.5 a 2.0 (ajuste conforme necess�rio)
    lb = [lambda_val, 0, 90, 90, 0.5, 0, 0];
    ub = [lambda_val, 0, 180, 180, 2.0, 3.0, 3.0];

    [x, fval, exitflag] = fmincon(@(X) fob(X, SVC, V, V_ini, alpha1, alpha2), ...
                                  x0, [], [], [], [], lb, ub, ...
                                  @(X) in(X, V, THDv, lii, SVC), options);

    if exitflag > 0
        lambda = x(1);
        lbd = x(2);
        alpha1 = x(3);
        alpha2 = x(4);
        FAKTOR = x(5:7);                 % +++++ atualiza fatores
        Ton1s = alpha1;                  % +++++ atualiza �ngulos base
        Ton2s = alpha2;
    else
        fprintf('Aviso: falha em lambda=%.3f (exitflag=%d). Mantendo valores anteriores.\n', lambda_val, exitflag);
    end

    Header_1;   % Atualiza V e THDv com os novos �ngulos e fatores

    % Armazena resultados
    LLLA(count) = lambda;
    % lbd_values(count) = lbd;
    ALPHA_1(count) = alpha1;
    ALPHA_2(count) = alpha2;
    FATORES(count, :) = FAKTOR;        % +++++ armazena fatores
    lambdas_armazenados(count) = lambda;

    for barra = 1:n_barras
        tensoes_3d(barra, :, count) = [V(barra, 3), V(barra, 4), V(barra, 5)];
        thdv_3d(barra, :, count)   = [THDv(barra, 1), THDv(barra, 2), THDv(barra, 3)];
    end

    lambda_anterior = lambda_val;
end
close(h_waitbar);

% ========== 6. SALVAR DADOS ==========
save('dados_analise.mat', 'LLLA', 'ALPHA_1', 'ALPHA_2', 'FATORES', 'lbd_values', ...
     'tensoes_3d', 'thdv_3d', 'lambdas_armazenados', 'n_barras');

fprintf('Dados salvos em dados_analise.mat\n');
fprintf('Tempo total: %.2f minutos\n', toc(t_ini)/60);

% ========== 7. PLOTAR RESULTADOS ==========
plotar_graficos_selecionados();