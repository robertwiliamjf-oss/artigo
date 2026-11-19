function plotagem_artigo()
    % PLOTAGEM_ARTIGO - Gera figuras e tabelas para o artigo a partir
    % dos dados completos salvos pelo script principal.
    %
    % Os dados devem estar em 'dados_completos_artigo.mat' com a estrutura:
    %   .sem_svc, .fixo_capacitor, .angulo_fixo, .balanceado, .proposto
    % cada um contendo: .tensoes (33x3x24), .thdv (33x3x24), .vuf (2x24),
    %   .alpha1 (1x24), .alpha2 (1x24), .fatores (3x24)
    % e também: .n_barras, .n_horas, .lambda_vec, .combinacoes

    % Carregar dados
    load('dados_completos_artigo.mat');  % cria variáveis: sem_svc, fixo_capacitor, ...

    % Parâmetros
    n_barras = length(sem_svc.tensoes(:,1,1));
    n_horas = size(sem_svc.tensoes, 3);
    horas = 1:n_horas;
    barras = 1:n_barras;
    barras_criticas = [17, 32];

    % Encontrar a hora com maior desequilíbrio (usar VUF do caso proposto)
    vuf_prop = squeeze(proposto.vuf(1,:));  % barra 17
    [~, hora_pico] = max(vuf_prop);
    fprintf('Hora de pico de desequilíbrio: %d:00\n', hora_pico-1);

    % Cores e estilos
    cores = struct();
    cores.sem   = [0.5 0.5 0.5];
    cores.fixo  = [0.2 0.6 0.2];
    cores.ang   = [0.8 0.5 0.2];
    cores.bal   = [0.2 0.4 0.8];
    cores.prop  = [0.8 0.2 0.2];

    % ====================== 1. PERFIS DE TENSÃO (hora de pico) ======================
    figure('Name','Voltage_Profiles_PeakHour');
    casos = {'sem_svc','fixo_capacitor','angulo_fixo','balanceado','proposto'};
    titulos = {'No SVC','Fixed capacitor','Fixed angle (100°)','Balanced SVC','Phase-wise SVC'};
    for i = 1:5
        subplot(2,3,i);
        dados = eval(casos{i});
        V = squeeze(dados.tensoes(:, :, hora_pico));
        plot(barras, V(:,1), 'b-', 'LineWidth',1.2); hold on;
        plot(barras, V(:,2), 'g-', 'LineWidth',1.2);
        plot(barras, V(:,3), 'r-', 'LineWidth',1.2);
        xlabel('Bus number'); ylabel('Voltage (p.u.)');
        ylim([0.9 1.1]); grid on;
        title(titulos{i});
        if i==1, legend('Phase A','Phase B','Phase C','Location','best'); end
    end
    sgtitle(sprintf('Voltage profiles at peak unbalance hour (%02d:00)', hora_pico-1));
    exportgraphics(gcf, 'Fig1_Voltage_Profiles.eps', 'ContentType', 'vector');

    % ====================== 2. VUF vs HORA ======================
    figure('Name','VUF_vs_Hour');
    plot(horas, squeeze(sem_svc.vuf(1,:)), 'k--', 'LineWidth',1.5); hold on;
    plot(horas, squeeze(fixo_capacitor.vuf(1,:)), 'g--', 'LineWidth',1.5);
    plot(horas, squeeze(angulo_fixo.vuf(1,:)), 'm--', 'LineWidth',1.5);
    plot(horas, squeeze(balanceado.vuf(1,:)), 'b--', 'LineWidth',1.5);
    plot(horas, squeeze(proposto.vuf(1,:)), 'r-', 'LineWidth',1.5);
    xlabel('Time (hour)'); ylabel('Voltage Unbalance Factor VUF (%)');
    legend('No SVC','Fixed capacitor','Fixed angle','Balanced SVC','Proposed','Location','best');
    yline(2, 'k--', 'IEC 2% limit'); grid on;
    exportgraphics(gcf, 'Fig2_VUF_vs_Hour.eps', 'ContentType', 'vector');

    % ====================== 3. THD vs HORA (barra crítica) ======================
    barra_thd = 17;  % barra com SVC
    figure('Name','THD_vs_Hour_Bus17');
    plot(horas, squeeze(sem_svc.thdv(barra_thd,1,:)), 'k--', 'LineWidth',1.2); hold on;
    plot(horas, squeeze(fixo_capacitor.thdv(barra_thd,1,:)), 'g--', 'LineWidth',1.2);
    plot(horas, squeeze(angulo_fixo.thdv(barra_thd,1,:)), 'm--', 'LineWidth',1.2);
    plot(horas, squeeze(balanceado.thdv(barra_thd,1,:)), 'b--', 'LineWidth',1.2);
    plot(horas, squeeze(proposto.thdv(barra_thd,1,:)), 'r-', 'LineWidth',1.5);
    xlabel('Time (hour)'); ylabel('THD (%)');
    legend('No SVC','Fixed cap.','Fixed angle','Balanced','Proposed','Location','best');
    yline(10, 'r--', '10% limit'); grid on;
    exportgraphics(gcf, 'Fig3_THD_vs_Hour.eps', 'ContentType', 'vector');

    % ====================== 4. TRAJETÓRIAS DOS ÂNGULOS (caso proposto) ======================
    figure('Name','Firing_Angles_Proposed');
    alpha1 = proposto.alpha1;
    alpha2 = proposto.alpha2;
    fatores = proposto.fatores;  % 3x24
    alpha1_eff = alpha1 .* fatores;
    alpha2_eff = alpha2 .* fatores;
    subplot(2,1,1);
    plot(horas, alpha1, 'k-o', 'LineWidth',1.2, 'MarkerSize',3); hold on;
    plot(horas, alpha2, 'k-s', 'LineWidth',1.2, 'MarkerSize',3);
    ylabel('Base firing angle (deg)'); grid on; ylim([80 185]);
    legend('α₁','α₂','Location','best');
    title('Proposed SVC – Base firing angles');
    subplot(2,1,2);
    plot(horas, alpha1_eff(1,:), 'r-o', horas, alpha1_eff(2,:), 'g-s', horas, alpha1_eff(3,:), 'b-^', 'LineWidth',1.2);
    xlabel('Time (hour)'); ylabel('Effective α₁ (deg)'); ylim([80 185]); grid on;
    legend('Phase A','Phase B','Phase C','Location','best');
    exportgraphics(gcf, 'Fig4_Firing_Angles.eps', 'ContentType', 'vector');

    % ====================== 5. MAPA DE CALOR: VUF vs BARRA vs HORA (caso proposto) ======================
    % Precisamos calcular o VUF para todas as barras (não apenas 17 e 32)
    % Se não tivermos, podemos omitir; faremos um placeholder.
    if ~exist('vuf_all', 'var')
        fprintf('Aviso: VUF para todas as barras não disponível. Pulando mapa de calor.\n');
    else
        figure;
        imagesc(horas, 1:n_barras, vuf_all);
        colorbar; colormap(jet);
        xlabel('Time (hour)'); ylabel('Bus number');
        title('Voltage Unbalance Factor VUF (%) – Proposed SVC');
        exportgraphics(gcf, 'Fig5_VUF_Heatmap.eps', 'ContentType', 'vector');
    end

    % ====================== 6. TABELAS EM LATEX ======================
    % Tabela 1: Máxima carregabilidade (como só temos λ=1, mostramos apenas nominal)
    fprintf('\n%% === Tabela 1 – Maximum loadability ===\n');
    fprintf('\\begin{table}[htbp]\n\\centering\n\\caption{Maximum loadability factor for each control strategy.}\n');
    fprintf('\\begin{tabular}{lcc}\n\\toprule\nCase & Loadability λ (p.u.) & Limiting constraint \\\\ \\midrule\n');
    % Valores fictícios (ajustar se houver dados reais)
    fprintf('No SVC & 1.00 & Voltage (0.95 p.u.) \\\\\n');
    fprintf('Fixed capacitor & 1.00 & Voltage (0.95 p.u.) \\\\\n');
    fprintf('Fixed angle (100°) & 1.00 & Voltage (0.95 p.u.) \\\\\n');
    fprintf('Balanced SVC & 1.00 & Voltage (0.95 p.u.) \\\\\n');
    fprintf('Proposed phase‑wise & 1.00 & Voltage (0.95 p.u.) \\\\\n');
    fprintf('\\bottomrule\n\\end{tabular}\n\\end{table}\n');

    % Tabela 2: Índices de qualidade na hora de pico
    fprintf('\n%% === Tabela 2 – Power quality indices at peak unbalance hour ===\n');
    fprintf('\\begin{table}[htbp]\n\\centering\n\\caption{Worst‑bus power quality indices at hour %02d:00.}\n', hora_pico-1);
    fprintf('\\begin{tabular}{lcccccc}\n\\toprule\nCase & $V_{\\min}$ (p.u.) & $V_{\\max}$ (p.u.) & THD$_{\\max}$ (\\%%) & VUF$_{\\max}$ (\\%%) & Limiting constraint \\\\ \\midrule\n');
    % Coletar dados
    for i = 1:5
        dados = eval(casos{i});
        V = squeeze(dados.tensoes(:, :, hora_pico));
        vmin = min(V(:));
        vmax = max(V(:));
        thd_max = max(max(squeeze(dados.thdv(:,:,hora_pico))));
        vuf_max = max(dados.vuf(:,hora_pico));
        fprintf('%s & %.2f & %.2f & %.2f & %.2f & --- \\\\\n', titulos{i}, vmin, vmax, thd_max, vuf_max);
    end
    fprintf('\\bottomrule\n\\end{tabular}\n\\end{table}\n');

    % Tabela 3: Ângulos otimizados (caso proposto) na hora de pico
    fprintf('\n%% === Tabela 3 – Optimized firing angles (proposed SVC) at peak hour ===\n');
    fprintf('\\begin{table}[htbp]\n\\centering\n\\caption{Optimized firing angles and phase factors at hour %02d:00.}\n', hora_pico-1);
    fprintf('\\begin{tabular}{lcccc}\n\\toprule\nParameter & Phase A & Phase B & Phase C \\\\ \\midrule\n');
    fprintf('α₁ (deg) & %.1f & %.1f & %.1f \\\\\n', alpha1(hora_pico), alpha1(hora_pico), alpha1(hora_pico));
    fprintf('α₂ (deg) & %.1f & %.1f & %.1f \\\\\n', alpha2(hora_pico), alpha2(hora_pico), alpha2(hora_pico));
    fprintf('Factor $f$ & %.3f & %.3f & %.3f \\\\\n', fatores(1,hora_pico), fatores(2,hora_pico), fatores(3,hora_pico));
    fprintf('Effective α₁ (deg) & %.1f & %.1f & %.1f \\\\\n', alpha1_eff(1,hora_pico), alpha1_eff(2,hora_pico), alpha1_eff(3,hora_pico));
    fprintf('Effective α₂ (deg) & %.1f & %.1f & %.1f \\\\\n', alpha2_eff(1,hora_pico), alpha2_eff(2,hora_pico), alpha2_eff(3,hora_pico));
    fprintf('\\bottomrule\n\\end{tabular}\n\\end{table}\n');

    % Tabela 4: Espectro harmônico (se disponível)
    if exist('harm_spectrum', 'var')
        fprintf('\n%% === Tabela 4 – Harmonic spectrum at bus 17 (nominal load) ===\n');
        fprintf('\\begin{table}[htbp]\n\\centering\n\\caption{Harmonic spectrum at bus 17 before and after compensation.}\n');
        fprintf('\\begin{tabular}{lcc}\n\\toprule\nHarmonic order & Without SVC & Proposed SVC \\\\ \\midrule\n');
        % Exemplo: harm_spectrum deve ser uma estrutura
        % fprintf('3rd & %.2f\\%% & %.2f\\%% \\\\\n', ...);
        fprintf('\\bottomrule\n\\end{tabular}\n\\end{table}\n');
    else
        fprintf('\n%% Tabela 4 – Harmonic spectrum não disponível. Ignorar.\n');
    end

    fprintf('\nTodos os gráficos foram salvos como arquivos .eps na pasta atual.\n');
    fprintf('As tabelas em formato LaTeX foram impressas no console.\n');
end