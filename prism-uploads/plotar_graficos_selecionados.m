function plotar_graficos_selecionados()
    
    load('dados_analise.mat', 'LLLA', 'ALPHA_1', 'ALPHA_2', ...
         'tensoes_3d', 'thdv_3d', 'lambdas_armazenados', 'n_barras');

    fprintf('\n--- MENU DE GRAFICOS ---\n');
    fprintf('1 - Tensao vs iteracao (curva diaria)\n');
    fprintf('2 - THDv vs iteracao (curva diaria)\n');
    fprintf('3 - Tensao vs lambda (original)\n');
    fprintf('4 - THDv vs lambda (original)\n');
    fprintf('5 - Perfil de tensao trifasico (hold on)\n');
    fprintf('6 - Perfil de THDv trifasico (hold on)\n');
    fprintf('7 - Angulos alpha1 e alpha2 vs iteracao\n');
    fprintf('8 - Todos os graficos acima\n');
    escolha = input('Digite o numero da opcao desejada: ');

    % Se escolheu opções 1,2,8, pedir barras
    if any(escolha == [1,2,8])
        fprintf('\nDigite os numeros das barras desejadas (ex: 1,5,10,17) ou "all" para todas: ');
        str_barras = input('', 's');
        if strcmpi(str_barras, 'all')
            barras_selecionadas = 1:n_barras;
        else
            barras_selecionadas = str2num(str_barras); %#ok<ST2NM>
            if isempty(barras_selecionadas)
                error('Entrada invalida. Use numeros separados por virgula ou "all".');
            end
        end
    end

    % 1. Tensão vs iteração (sequencial)
    if escolha == 1 || escolha == 8
        figura1 = figure('Name','Tensao vs Iteracao (curva diaria)');
        hold on;
        cores = lines(length(barras_selecionadas));
        for idx = 1:length(barras_selecionadas)
            barra = barras_selecionadas(idx);
            tensao_media = squeeze(mean(tensoes_3d(barra,:,:), 2));
            plot(1:length(tensao_media), tensao_media, ...
                 'Color', cores(idx,:), 'LineWidth', 1.5, ...
                 'DisplayName', sprintf('Barra %d', barra));
        end
        xlabel('Iteracao (ordem de carregamento)');
        ylabel('Tensao (p.u.)');
        title('Tensao ao longo da sequencia de carregamento');
        legend show; grid on;

        % Adicionar rótulos dos lambdas no eixo x (opcional)
        xticks(1:length(lambdas_armazenados));
        xticklabels(arrayfun(@(x) sprintf('%.2f', x), lambdas_armazenados, 'UniformOutput', false));
        xlabel('Iteracao (valores de lambda indicados abaixo)');
        % Girar rótulos se necessário
        xtickangle(45);
    end

    % 2. THDv vs iteração
    if escolha == 2 || escolha == 8
        figura2 = figure('Name','THDv vs Iteracao (curva diaria)');
        hold on;
        cores = lines(length(barras_selecionadas));
        for idx = 1:length(barras_selecionadas)
            barra = barras_selecionadas(idx);
            thdv_media = squeeze(mean(thdv_3d(barra,:,:), 2));
            plot(1:length(thdv_media), thdv_media, ...
                 'Color', cores(idx,:), 'LineWidth', 1.5, ...
                 'DisplayName', sprintf('Barra %d', barra));
        end
        xlabel('Iteracao (ordem de carregamento)');
        ylabel('THDv (%)');
        title('THDv ao longo da sequencia de carregamento');
        legend show; grid on;

        xticks(1:length(lambdas_armazenados));
        xticklabels(arrayfun(@(x) sprintf('%.2f', x), lambdas_armazenados, 'UniformOutput', false));
        xtickangle(45);
    end

    % 3. Tensão vs lambda (original)
    if escolha == 3 || escolha == 8
        % ... código anterior ...
    end

    % 4. THDv vs lambda (original)
    if escolha == 4 || escolha == 8
        % ... código anterior ...
    end

    % 5. Perfil de tensão trifásico (hold on)
    if escolha == 5 || escolha == 8
        % ... código anterior ...
    end

    % 6. Perfil de THDv trifásico (hold on)
    if escolha == 6 || escolha == 8
        % ... código anterior ...
    end

    % 7. Ângulos vs iteração
    if escolha == 7 || escolha == 8
        figura7 = figure('Name','Angulos vs Iteracao');
        plot(1:length(ALPHA_1), ALPHA_1, 'ro-', 'LineWidth', 1.5, 'DisplayName', '\alpha_1');
        hold on;
        plot(1:length(ALPHA_2), ALPHA_2, 'bs-', 'LineWidth', 1.5, 'DisplayName', '\alpha_2');
        xlabel('Iteracao');
        ylabel('Angulo (graus)');
        title('Angulos de disparo ao longo da sequencia de carregamento');
        legend show; grid on;

        xticks(1:length(lambdas_armazenados));
        xticklabels(arrayfun(@(x) sprintf('%.2f', x), lambdas_armazenados, 'UniformOutput', false));
        xtickangle(45);
    end

    fprintf('\nGraficos gerados.\n');
end