function plotar_resultados_artigo()
clc
close all
load('dados_analise.mat')
hora = 1:24;
barra = 17;
%% ==========================================================
% FIGURE 1 - DAILY LOADING PROFILES
%% ==========================================================
figure('Name','Daily Loading Profiles');
hold on
plot(hora,des1_vec,'LineWidth',2)
plot(hora,des2_vec,'LineWidth',2)
plot(hora,des3_vec,'LineWidth',2)
xlabel('Hour')
ylabel('Loading factor')
legend('Phase A','Phase B','Phase C',...
       'Location','best')
ytickformat('%.4f') % Força 4 casas decimais no eixo Y
grid on
%% ==========================================================
% FIGURE 2 - VUF
%% ==========================================================
VUF_opt = squeeze(VUF_hist(1,1,:));
figure('Name','Voltage Unbalance Factor');
hold on
plot(hora,VUF_fixo(1,:),'--','LineWidth',2,...
    'DisplayName','Fixed SVC')
plot(hora,VUF_opt,'LineWidth',2,...
    'DisplayName','Optimized SVC')
xlabel('Hour')
ylabel('VUF (%)')
legend('show','Location','best')
ytickformat('%.4f') % Força 4 casas decimais no eixo Y
grid on
%% ==========================================================
% FIGURE 3 - NEGATIVE SEQUENCE VOLTAGE
%% ==========================================================
V2_fix = zeros(1,24);
V2_opt = zeros(1,24);
a = exp(1j*2*pi/3);
for h = 1:24
    %% FIXED
    Va = tensoes_fixo(barra,1,h) ...
        * exp(1j*T_fixo(barra,1,h));
    Vb = tensoes_fixo(barra,2,h) ...
        * exp(1j*T_fixo(barra,2,h));
    Vc = tensoes_fixo(barra,3,h) ...
        * exp(1j*T_fixo(barra,3,h));
    V2 = (Va + a^2*Vb + a*Vc)/3;
    V2_fix(h) = abs(V2);
    %% OPTIMIZED
    Va = tensoes_4d(barra,1,1,h) ...
        * exp(1j*T_4d(barra,1,1,h));
    Vb = tensoes_4d(barra,2,1,h) ...
        * exp(1j*T_4d(barra,2,1,h));
    Vc = tensoes_4d(barra,3,1,h) ...
        * exp(1j*T_4d(barra,3,1,h));
    V2 = (Va + a^2*Vb + a*Vc)/3;
    V2_opt(h) = abs(V2);
end
figure('Name','Negative Sequence Voltage')
hold on
plot(hora,V2_fix,'--','LineWidth',2,...
    'DisplayName','Fixed SVC')
plot(hora,V2_opt,'LineWidth',2,...
    'DisplayName','Optimized SVC')
xlabel('Hour')
ylabel('|V_2| (p.u.)')
legend('show','Location','best')
ytickformat('%.4f') % Força 4 casas decimais no eixo Y
grid on
%% ==========================================================
% FIGURE 4 - BUS 32 VOLTAGES
%% ==========================================================
figure('Name','Bus Voltages')
hold on
plot(hora,squeeze(tensoes_fixo(barra,1,:)),...
    '--','LineWidth',2,...
    'DisplayName','Phase A - Fixed')
plot(hora,squeeze(tensoes_4d(barra,1,1,:)),...
    'LineWidth',2,...
    'DisplayName','Phase A - Optimized')
plot(hora,squeeze(tensoes_fixo(barra,2,:)),...
    '--','LineWidth',2,...
    'DisplayName','Phase B - Fixed')
plot(hora,squeeze(tensoes_4d(barra,2,1,:)),...
    'LineWidth',2,...
    'DisplayName','Phase B - Optimized')
plot(hora,squeeze(tensoes_fixo(barra,3,:)),...
    '--','LineWidth',2,...
    'DisplayName','Phase C - Fixed')
plot(hora,squeeze(tensoes_4d(barra,3,1,:)),...
    'LineWidth',2,...
    'DisplayName','Phase C - Optimized')
xlabel('Hour')
ylabel('Voltage magnitude (p.u.)')
legend('show','Location','best')
ylim([0.9 1.1])
ytickformat('%.4f') % Força 4 casas decimais no eixo Y
grid on
%% ==========================================================
% FIGURE 5 - EFFECTIVE FIRING ANGLES
%% ==========================================================
alphaA = squeeze(ALPHA_1(1,:)) ...
       .* squeeze(FATORES(1,1,:))';
alphaB = squeeze(ALPHA_1(1,:)) ...
       .* squeeze(FATORES(1,2,:))';
alphaC = squeeze(ALPHA_1(1,:)) ...
       .* squeeze(FATORES(1,3,:))';
figure('Name','Effective Firing Angles')
hold on
plot(hora,alphaA,'b-o', 'LineWidth', 1.5, 'MarkerSize', 4,...
    'DisplayName','Phase A')
plot(hora,alphaB,'g-o', 'LineWidth', 1.5, 'MarkerSize', 4,...
    'DisplayName','Phase B')
plot(hora,alphaC,'r-o', 'LineWidth', 1.5, 'MarkerSize', 4,...
    'DisplayName','Phase C')
xlabel('Hour')
ylabel('Effective firing angle (deg)')
legend('show','Location','best')
ytickformat('%.4f') % Força 4 casas decimais no eixo Y
grid on
%% ==========================================================
% FIGURE 6 - PHASORS AT WORST CONDITION
%% ==========================================================
[~,hora_critica] = max(VUF_fixo(1,:));
Va_fix = tensoes_fixo(barra,1,hora_critica) ...
    * exp(1j*T_fixo(barra,1,hora_critica));
Vb_fix = tensoes_fixo(barra,2,hora_critica) ...
    * exp(1j*T_fixo(barra,2,hora_critica));
Vc_fix = tensoes_fixo(barra,3,hora_critica) ...
    * exp(1j*T_fixo(barra,3,hora_critica));
Va_opt = tensoes_4d(barra,1,1,hora_critica) ...
    * exp(1j*T_4d(barra,1,1,hora_critica));
Vb_opt = tensoes_4d(barra,2,1,hora_critica) ...
    * exp(1j*T_4d(barra,2,1,hora_critica));
Vc_opt = tensoes_4d(barra,3,1,hora_critica) ...
    * exp(1j*T_4d(barra,3,1,hora_critica));
figure('Name','Voltage Phasors')
hold on
axis equal
grid on
quiver(0,0,real(Va_fix),imag(Va_fix),0,...
    '--','LineWidth',2)
quiver(0,0,real(Vb_fix),imag(Vb_fix),0,...
    '--','LineWidth',2)
quiver(0,0,real(Vc_fix),imag(Vc_fix),0,...
    '--','LineWidth',2)
quiver(0,0,real(Va_opt),imag(Va_opt),0,...
    'LineWidth',2)
quiver(0,0,real(Vb_opt),imag(Vb_opt),0,...
    'LineWidth',2)
quiver(0,0,real(Vc_opt),imag(Vc_opt),0,...
    'LineWidth',2)
xlabel('Real axis')
ylabel('Imaginary axis')
legend('Phase A - Fixed',...
       'Phase B - Fixed',...
       'Phase C - Fixed',...
       'Phase A - Optimized',...
       'Phase B - Optimized',...
       'Phase C - Optimized',...
       'Location','best')
xtickformat('%.4f') % Força 4 casas decimais no eixo X (Fasores)
ytickformat('%.4f') % Força 4 casas decimais no eixo Y (Fasores)
fprintf('\nWorst operating hour = %d\n',hora_critica);
%% ==========================================================
% FIGURE 7 - VOLTAGE DEVIATION FROM PHASE AVERAGE
%% ==========================================================
Va_fix = squeeze(tensoes_fixo(barra,1,:));
Vb_fix = squeeze(tensoes_fixo(barra,2,:));
Vc_fix = squeeze(tensoes_fixo(barra,3,:));
Vm_fix = (Va_fix + Vb_fix + Vc_fix)/3;
dA_fix = Va_fix - Vm_fix;
dB_fix = Vb_fix - Vm_fix;
dC_fix = Vc_fix - Vm_fix;
Va_opt = squeeze(tensoes_4d(barra,1,1,:));
Vb_opt = squeeze(tensoes_4d(barra,2,1,:));
Vc_opt = squeeze(tensoes_4d(barra,3,1,:));
Vm_opt = (Va_opt + Vb_opt + Vc_opt)/3;
dA_opt = Va_opt - Vm_opt;
dB_opt = Vb_opt - Vm_opt;
dC_opt = Vc_opt - Vm_opt;
figure('Name','Voltage Deviation From Average')
subplot(3,1,1)
plot(hora,dA_fix,'--','LineWidth',2,...
    'DisplayName','Fixed SVC')
hold on
plot(hora,dA_opt,'LineWidth',2,...
    'DisplayName','Optimized SVC')
ylabel('\DeltaV_A (p.u.)')
legend('show','Location','best')
ytickformat('%.4f')
grid on
subplot(3,1,2)
plot(hora,dB_fix,'--','LineWidth',2,...
    'DisplayName','Fixed SVC')
hold on
plot(hora,dB_opt,'LineWidth',2,...
    'DisplayName','Optimized SVC')
ylabel('\DeltaV_B (p.u.)')
legend('show','Location','best')
ytickformat('%.4f')
grid on
subplot(3,1,3)
plot(hora,dC_fix,'--','LineWidth',2,...
    'DisplayName','Fixed SVC')
hold on
plot(hora,dC_opt,'LineWidth',2,...
    'DisplayName','Optimized SVC')
ylabel('\DeltaV_C (p.u.)')
xlabel('Hour')
legend('show','Location','best')
ytickformat('%.4f')
grid on
%% ==========================================================
% FIGURE 8 - VOLTAGE DEVIATION INDEX
%% ==========================================================
DV_fix = sqrt( ...
    dA_fix.^2 + ...
    dB_fix.^2 + ...
    dC_fix.^2 );
DV_opt = sqrt( ...
    dA_opt.^2 + ...
    dB_opt.^2 + ...
    dC_opt.^2 );
figure('Name','Voltage Deviation Index')
plot(hora,DV_fix,'--','LineWidth',2,...
    'DisplayName','Fixed SVC')
hold on
plot(hora,DV_opt,'LineWidth',2,...
    'DisplayName','Optimized SVC')
xlabel('Hour')
ylabel('Voltage deviation index (p.u.)')
legend('show','Location','best')
ytickformat('%.4f')
grid on
%% ==========================================================
% FIGURE 9 - VOLTAGE ENVELOPE
%% ==========================================================
barra = 32;
hora = 1:24;
%% FIXED CASE
Va_fix = squeeze(tensoes_fixo(barra,1,:));
Vb_fix = squeeze(tensoes_fixo(barra,2,:));
Vc_fix = squeeze(tensoes_fixo(barra,3,:));
Vmax_fix = max([Va_fix Vb_fix Vc_fix],[],2);
Vmin_fix = min([Va_fix Vb_fix Vc_fix],[],2);
Vm_fix   = (Va_fix + Vb_fix + Vc_fix)/3;
%% OPTIMIZED CASE
Va_opt = squeeze(tensoes_4d(barra,1,1,:));
Vb_opt = squeeze(tensoes_4d(barra,2,1,:));
Vc_opt = squeeze(tensoes_4d(barra,3,1,:));
Vmax_opt = max([Va_opt Vb_opt Vc_opt],[],2);
Vmin_opt = min([Va_opt Vb_opt Vc_opt],[],2);
Vm_opt   = (Va_opt + Vb_opt + Vc_opt)/3;
%% PLOT
figure('Name','Voltage Envelope');
hold on
box on
% Envelope - Fixed
fill([hora fliplr(hora)], ...
     [Vmax_fix' fliplr(Vmin_fix')], ...
     [0.85 0.85 0.85], ...
     'EdgeColor','none', ...
     'FaceAlpha',0.60, ...
     'DisplayName','Fixed SVC Envelope');
% Envelope - Optimized
fill([hora fliplr(hora)], ...
     [Vmax_opt' fliplr(Vmin_opt')], ...
     [0.30 0.60 1.00], ...
     'EdgeColor','none', ...
     'FaceAlpha',0.35, ...
     'DisplayName','Optimized SVC Envelope');
% Average voltage
plot(hora,Vm_fix,...
    'k--',...
    'LineWidth',2,...
    'DisplayName','Average Voltage - Fixed');
plot(hora,Vm_opt,...
    'b',...
    'LineWidth',2,...
    'DisplayName','Average Voltage - Optimized');
% Upper/lower limits
plot(hora,Vmax_fix,...
    'k:',...
    'LineWidth',1);
plot(hora,Vmin_fix,...
    'k:',...
    'LineWidth',1);
plot(hora,Vmax_opt,...
    'b:',...
    'LineWidth',1);
plot(hora,Vmin_opt,...
    'b:',...
    'LineWidth',1);
xlabel('Hour')
ylabel('Voltage magnitude (p.u.)')
legend('Location','best')
ytickformat('%.4f')
grid on
% Ajuste automático da escala vertical
margem = 0.005;
ylim([ ...
    min([Vmin_fix;Vmin_opt]) - margem , ...
    max([Vmax_fix;Vmax_opt]) + margem ...
]);
xlim([1 24])
%% ==========================================================
% FIGURE 10 - THD ENVELOPE
%% ==========================================================
barra = 32;
hora = 1:24;
%% FIXED CASE
THDa_fix = squeeze(thdv_fixo(barra,1,:));
THDb_fix = squeeze(thdv_fixo(barra,2,:));
THDc_fix = squeeze(thdv_fixo(barra,3,:));
THDmax_fix = max([THDa_fix THDb_fix THDc_fix],[],2);
THDmin_fix = min([THDa_fix THDb_fix THDc_fix],[],2);
THDavg_fix = (THDa_fix + THDb_fix + THDc_fix)/3;
%% OPTIMIZED CASE
THDa_opt = squeeze(thdv_4d(barra,1,1,:));
THDb_opt = squeeze(thdv_4d(barra,2,1,:));
THDc_opt = squeeze(thdv_4d(barra,3,1,:));
THDmax_opt = max([THDa_opt THDb_opt THDc_opt],[],2);
THDmin_opt = min([THDa_opt THDb_opt THDc_opt],[],2);
THDavg_opt = (THDa_opt + THDb_opt + THDc_opt)/3;
%% PLOT
figure('Name','THD Envelope')
hold on
box on
% Envelope - Fixed SVC
fill([hora fliplr(hora)], ...
     [THDmax_fix' fliplr(THDmin_fix')], ...
     [0.85 0.85 0.85], ...
     'EdgeColor','none', ...
     'FaceAlpha',0.60, ...
     'DisplayName','Fixed SVC Envelope');
% Envelope - Optimized SVC
fill([hora fliplr(hora)], ...
     [THDmax_opt' fliplr(THDmin_opt')], ...
     [0.30 0.60 1.00], ...
     'EdgeColor','none', ...
     'FaceAlpha',0.35, ...
     'DisplayName','Optimized SVC Envelope');
% Average THD
plot(hora,THDavg_fix,...
    'k--',...
    'LineWidth',2,...
    'DisplayName','Average THD - Fixed');
plot(hora,THDavg_opt,...
    'b',...
    'LineWidth',2,...
    'DisplayName','Average THD - Optimized');
% Upper/Lower bounds
plot(hora,THDmax_fix,...
    'k:',...
    'LineWidth',1);
plot(hora,THDmin_fix,...
    'k:',...
    'LineWidth',1);
plot(hora,THDmax_opt,...
    'b:',...
    'LineWidth',1);
plot(hora,THDmin_opt,...
    'b:',...
    'LineWidth',1);
% PRODIST limit (ajuste se desejar outro valor)
yline(8,...
      'r--',...
      'LineWidth',1.5,...
      'DisplayName','PRODIST Limit');
xlabel('Hour')
ylabel('THD_V (%)')
legend('Location','best')
ytickformat('%.4f')
grid on
xlim([1 24])
% Ajuste automático do eixo Y
margem = 0.5;
ylim([ ...
    max(0,min([THDmin_fix;THDmin_opt])-margem), ...
    max([THDmax_fix;THDmax_opt])+margem ...
]);
%% ==========================================================
% FIGURE 11 - THD PER PHASE
%% ==========================================================
barra = 32;
hora = 1:24;
THDa = squeeze(thdv_4d(barra,1,1,:));
THDb = squeeze(thdv_4d(barra,2,1,:));
THDc = squeeze(thdv_4d(barra,3,1,:));
figure('Name','THD per Phase')
plot(hora,THDa,...
    'LineWidth',2,...
    'DisplayName','Phase A')
hold on
plot(hora,THDb,...
    'LineWidth',2,...
    'DisplayName','Phase B')
plot(hora,THDc,...
    'LineWidth',2,...
    'DisplayName','Phase C')
xlabel('Hour')
ylabel('THD_V (%)')
legend('show','Location','best')
ytickformat('%.4f')
grid on
xlim([1 24])
%% ==========================================================
% FIGURA 12 - TENSÃO DE SEQUÊNCIA ZERO
%% ==========================================================
V0_fix = zeros(1,24);
V0_opt = zeros(1,24);
for h = 1:24
    %% FIXED
    Va = tensoes_fixo(barra,1,h) ...
        * exp(1j*T_fixo(barra,1,h));
    Vb = tensoes_fixo(barra,2,h) ...
        * exp(1j*T_fixo(barra,2,h));
    Vc = tensoes_fixo(barra,3,h) ...
        * exp(1j*T_fixo(barra,3,h));
    V0 = (Va + Vb + Vc)/3;          % sequência zero
    V0_fix(h) = abs(V0);
    %% OPTIMIZED
    Va = tensoes_4d(barra,1,1,h) ...
        * exp(1j*T_4d(barra,1,1,h));
    Vb = tensoes_4d(barra,2,1,h) ...
        * exp(1j*T_4d(barra,2,1,h));
    Vc = tensoes_4d(barra,3,1,h) ...
        * exp(1j*T_4d(barra,3,1,h));
    V0 = (Va + Vb + Vc)/3;          % sequência zero
    V0_opt(h) = abs(V0);
end
figure('Name','Zero Sequence Voltage')
hold on
plot(hora, V0_fix, '--', 'LineWidth', 2, ...
    'DisplayName', 'Fixed SVC')
plot(hora, V0_opt, 'LineWidth', 2, ...
    'DisplayName', 'Optimized SVC')
xlabel('Hour')
ylabel('|V_0| (p.u.)')
legend('show', 'Location', 'best')
ytickformat('%.4f')
grid on
%% ==========================================================
% SUMMARY
%% ==========================================================
fprintf('\n');
fprintf('Maximum VUF (Fixed)      = %.4f %%\n',max(VUF_fixo(1,:)));
fprintf('Maximum VUF (Optimized)  = %.4f %%\n',max(VUF_opt));
fprintf('Average VUF (Fixed)      = %.4f %%\n',mean(VUF_fixo(1,:)));
fprintf('Average VUF (Optimized)  = %.4f %%\n',mean(VUF_opt));
fprintf('\nAll article figures generated successfully.\n');
fprintf('\n');
% Alterado de %.5f para %.4f para manter o padrão estrito de 4 casas decimais
fprintf('Maximum DV (Fixed)     = %.4f pu\n',max(DV_fix));
fprintf('Maximum DV (Optimized) = %.4f pu\n',max(DV_opt));
fprintf('Average DV (Fixed)     = %.4f pu\n',mean(DV_fix));
fprintf('Average DV (Optimized) = %.4f pu\n',mean(DV_opt));
end