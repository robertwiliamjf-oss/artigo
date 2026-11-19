% %        Unifica��o das Fun��es MICTH e HVDC (Em Constru��o)!       % 
% % % % Metodo de Inje�ao de Correntes Trifasico Harmonico (MICTH)% % %
% % % %  Aplicado `a Modelagem do Elo de Corrente Continua (HVDC) % % %
% %-------------------------------------------------------------------%
% 
% 
% % % EM PROCESSO DE TESTES. EM CONSTRUCAO!!! Acertar o Calculo de Toff!!!
% 
% % Como diferncial para a versao 'HVDC_10' est�o:
%       _A possibilidade de Analise de RCTs e SVCs (S�NCRONO APENAS!!!);
%       _Incorpora��o do calculo de cruzamentos por Zero da Corrente (EM CONSTRU��O E TESTES!);
%       _Cargas no modelo ZIP (CONFERIR!!!);
%       _Possibilidade de Utiliza��o de Fontes de Correntes Harm�nicas;
%       _Elimina��o das Barras VTeta em cada fun��o (Novas Vers�es), ficando
%   apenas a coloca��o de '1' na diagonal para ser feito nessa fun��o;
%       _Foi realizada tamb�m a substitui��o de 'Conv_6' para 'Conv_7', que
%   realiza a montagem da matriz Jacobiana de forma Vetorizada, aumentando
%   a Velocidade.
% 
% Os RCTs e SVCs est�o configurados para existirem apenas em SISTEMAS S�NCRONOS!!!
% (� preciso fazer da mesma forma que o feito para separar as barras e o
% mesmo que foi feito para os Conversores!!!)
% 
% 
% A inicializa��o do Vetor de Tens�es agora tb pode ser feito com 'Vret'
% contendo um n�mero diferente de harm�nicas armazenado, desde que o n�mero
% de Barras seja o mesmo. Na pr�tica isso permite que uma simula��o com um
% n�mero pequeno de harm�nicas possa ser utilizado como inicializa��o para
% simula��es com grande n�mero de harm�nicas, ou vice-versa.
% 
% A partir de agora a refer�ncia adotada � a FASORIAL, o que faz necess�rio
% um tratamento nos dados de entrada e sa�da:
%       Fourier  = [V/__+Fi = V*cos(nwt - Fi) = V*cos(Fi) + j*V*sen(Fi) = Vr + j*Vm, sendo Fi = atan2(Vm,Vr)].
%       Fasorial = [V/__+Fi = V*cos(nwt + Fi) = V*cos(Fi) + j*V*sen(Fi) = Vr + j*Vm, sendo Fi = atan2(Vm,Vr)].
%       Como EXEMPLO tomemos a funcao 'sen(nwt) = cos(nwt-90)':
%           Para Fourier:  'Vr = cos(+90) = 0' e 'Vm = sen(+90) = 1', ou seja, '0 + j*1' (Fi = 90), mas a referencia CORRETA seria:
%           Fasorialmente: 'Vr = cos(-90) = 0' e 'Vm = sen(-90) =-1',
%           ou seja, '0 - j*1' (Fi =-90). (Ocorre uma rota��o de 180� na parte Imagin�ria apenas)


%-------------------------------------------------------------------%
% % % Metodo de Inje�ao de Correntes Trifasico Harmonico (MICTH)% % %
% % %  % %
%-------------------------------------------------------------------%

% profile on
% clc

%%skip = 0                                    %%'1'   Entra com todas as condi��es iniciais (pode ser utilizado para retomar uma simulacao interrompida),
    % %                                             %%'0.5' Entra com apenas 'Vret' (para condi��es iniciais de Tens�o) e
    % %                                             %%'0'   Roda Normalmente
if skip==1
    cmax = cmax+iter;
    tempo1 = tempo;
    tempo=0;
    iter = iter+1;                          %%Por causa da localiza��o dos incrementos do contador
    tic
end
if skip==0 || skip==0.5

tempo=0;
tempo1=0;                                   %%Apenas para compatibilizar
tic

%% INICIALIZA�OES GERAIS
%%>
%%Feito em 'Header_x'
%%<

%% Carregando os Dados
%%>
%%Feito em 'Header_x' (Consultar 'HVDC_9' para os arquivos originalmente utilizados)
%%<
%% Leitura de Dados:
LerDados_3                                  %%Leitura de Dados
Num_Interna_2                               %%Reordenacao dos Numeros das Barras (Se preciso apenas)


iF = find(h==hf1);                          %%Numero da harmonica referente ao Bloco da Fundamental na Freq 1
F1 = 6*(iF-1)*Nb+1;                         %%Endere�o inicial da Fundamental na Freq 1

if Ent(1)==1
	iF = find(h==hf2);                      %%Numero da harmonica referente ao Bloco da Fundamental na Freq 2
	F2 = 6*(iF-1)*Nb+1;                     %%Endere�o inicial da Fundamental na Freq 2
end

%%Matriz para identifica��o das Harm�nicas e das Barras em fun��o da
%%Posi��o no vetor de Dados (apenas para an�lises externas):
%>
Barras = zeros(6*H*Nb,2);
z = 1:6*H*Nb;
Barras(:,1) = fix((z -(6*Nb*(fix(z/(6*Nb))+1-1) +1 -6))/6).';   %%Numero da Barra
z = (1:6*H*Nb)-1e-3;
Barras(:,2) = h(fix(z/(6*Nb))+1).';                             %%Freq Harm. que a Barra pertence

Barras(Barras(:,1)==0,1)=Nb;                                    %%Consertando para os Extremos
%<
% ======================================================================= %

% ======================================================================= %

% ======================================================================= %
%% INICIALIZA�AO DAS MATRIZES A SEREM UTILIZADAS
%%>
	xx1 = find(BarCtr~=0);                  %%Encontra as Barras que nao sao VTeta
    Lxx1= length(xx1);                      %%Numero de Barras com Tensao Controlada
    
    br = find(BarCtr==0);
    br = kg(br);                                %%N�mero das Barras de Referencia
    if Ent(1)==1
        xx2 = find(BarCtr2~=0);                 %%Encontra as Barras que nao sao VTeta no Sistema 2
        Lxx2= length(xx2);                      %%Numero de Barras com Tensao Controlada no Sistema 2
        Lxx = Lxx1 + Lxx2;
        
        br2 = find(BarCtr2==0);
        br2 = kg2(br2);                         %%N�mero das Barras de Referencia
        br = [br;br2];                          %%Todas as Barras de Refer�ncia
    else
        Lxx = Lxx1;
    end
    
    Ne = 3*(Lxx);                           %% Numero de eq. Adicionais (= numero de Geradores PV)
   
    % ------------------------------------
    J = sparse(6*H*Nb + Ne , 6*H*Nb + Ne);  %% Matriz Principal: J = -Jacobiano !!!
    J0= sparse(6*H*Nb + Ne , 6*H*Nb + Ne);  %% Matriz Inicial, Nao-Atualizavel!
    % ------------------------------------
    
    Dim = 6*H*Nb + Ne;
    %%fprintf('\nDimens�o Total do Sistema: %g x %g\n\n',Dim,Dim);
    
    
    %%Tive que trocar a posi��o do c�digo abaixo para o final para permitir
    %%que a Inicializa��o Pr�via tamb�m possa ser feita para cargas ZIP.
    %%Caso contr�rio as Matrizes 'Ycz' e 'Ycz_ref' seriam montadas de forma
    %%errada.
    %>
    % %     if skip==0                              %%Sem Inicializa��o Pr�via
    % %         Vret = zeros(6*H*Nb,1);
    % %     elseif skip==0.5                        %%Com Inicializa��o Pr�via (Total ou Parcial)
    % %         Vret = zeros(6*H*Nb,1);
    % %         if length(Vret)>=length(V)
    % %             Vret(1:length(V),1) = V;
    % %         else
    % %             Vret = V(1:length(Vret),1);
    % %         end
    % %     end
    
    Vret = zeros(6*H*Nb,1);
    %<

    DI  = ones(6*H*Nb + Ne,1);          	%%Vetor de Residuos de Corrente ('1' para iniciar o Loop Principal precocemente)
    DV  = zeros(Ne,1);                  	%%Residuos de Tensao das eq. de Controle (Modulo)
    if Ent(12)==1
        DInj  = zeros(6*H*Nb,1);           	%%Para as Fontes de Correntes Harm�nicas
    end
    
    Vdata = sparse(6*H*Nb,cmax);         	%%Armazenamento dos Dados
%%<

% % ========================= % % ========================= % % ========================= % %
% % ========================= % % ========================= % % ========================= % %

%% Calculando as Partes Reais e Imaginarias da Componente Fundamental da Tensao
%%>
if skip==0 || (Ent(5)==2 || Ent(5)==3)      %%(Trocado por causa das Cargas ZIP)
    %>
    if ABC==1                               %%Para Arquivo de Dados com Sequ�ncia de Fase ABC
    
    Var = Va.*cos(Ta);	Vam = Va.*sin(Ta);
    Vbr = Vb.*cos(Tb);	Vbm = Vb.*sin(Tb);
    Vcr = Vc.*cos(Tc);	Vcm = Vc.*sin(Tc);
    
    elseif ABC==0                           %%Para Arquivos de Dados com Sequ�ncia de Fase ACB (Arquivos Antigos sem YD)
    
    %%Dessa forma os �ngulos das Defasagens dos Transf. YD devem ser
    %%trocados em rela��o aos Arquivos de Dados Antigos (Tem que Incluir
    %%a Defasagem de -30� e n�o de +30)
    Var = Va.*cos(Ta);	Vam = Va.*sin(Ta);
    Vbr = Vb.*cos(Tc);	Vbm = Vb.*sin(Tc);
    Vcr = Vc.*cos(Tb);	Vcm = Vc.*sin(Tb);
    
    elseif ABC==-1                          %%Para Arquivos de Dados com Sequ�ncia de Fase ACB (Arquivos Antigos com YD)
    
    %%Dessa forma os �ngulos devem ser os mesmos dos Arquivos de Dados
    %%Antigos (Por causa das Defasagens dos Transf. YD)
    Var = Va.*cos(Ta);	Vam = Va.*sin(-Ta);
    Vbr = Vb.*cos(Tb);	Vbm = Vb.*sin(-Tb);
    Vcr = Vc.*cos(Tc);	Vcm = Vc.*sin(-Tc);
    
    end
    %<

    if Ent(1)==1
        %>
        if ABC==1                               %%Para Arquivo de Dados com Sequ�ncia de Fase ABC
        
        Var2 = Va2.*cos(Ta2);	Vam2 = Va2.*sin(Ta2);
        Vbr2 = Vb2.*cos(Tb2);	Vbm2 = Vb2.*sin(Tb2);
        Vcr2 = Vc2.*cos(Tc2);	Vcm2 = Vc2.*sin(Tc2);
        
        elseif ABC==0                           %%Para Arquivos de Dados com Sequ�ncia de Fase ACB (Arquivos Antigos sem YD)

            %%Dessa forma os �ngulos das Defasagens dos Transf. YD devem ser
            %%trocados em rela��o aos Arquivos de Dados Antigos (Tem que Incluir
            %%a Defasagem de -30� e n�o de +30)
            Var2 = Va2.*cos(Ta2);	Vam2 = Va2.*sin(Ta2);
            Vbr2 = Vb2.*cos(Tc2);	Vbm2 = Vb2.*sin(Tc2);
            Vcr2 = Vc2.*cos(Tb2);	Vcm2 = Vc2.*sin(Tb2);
            
        elseif ABC==-1                          %%Para Arquivos de Dados com Sequ�ncia de Fase ACB (Arquivos Antigos com YD)

            %%Dessa forma os �ngulos devem ser os mesmos dos Arquivos de Dados
            %%Antigos (Por causa das Defasagens dos Transf. YD)
            Var2 = Va2.*cos(Ta2);	Vam2 = Va2.*sin(-Ta2);
            Vbr2 = Vb2.*cos(Tb2);	Vbm2 = Vb2.*sin(-Tb2);
            Vcr2 = Vc2.*cos(Tc2);	Vcm2 = Vc2.*sin(-Tc2);            
        end
        %<
    end
%%<

%% Montagem dos Vetores de Tensao Para a Freq Fundamental
%%Cada bloco contem as componentes das 3 Fases
%%>
    if Ent(1)==0, Nb1 = Nb; end                      %%Apenas para Compatibilizar
    for k = 1:Nb1
        i = F1 + 6*kb(k) - 6;                   %%Cada bloco contem 6 linhas (Saltando as componentes anteriores)

        Vret(i:i+5,1) = [Var(k) ; Vbr(k) ; Vcr(k) ;     	%%Parte Real das Tensoes das 3 Fases
                         Vam(k) ; Vbm(k) ; Vcm(k)];         %%Parte Imaginaria das Tensoes das 3 Fases
    end
    if Ent(1)==1
        for k = 1:Nb2
            i = F2 + 6*kb2(k) - 6;              %%Cada bloco contem 6 linhas (Saltando as componentes anteriores)

            Vret(i:i+5,1) = [Var2(k) ; Vbr2(k) ; Vcr2(k) ;  %%Parte Real das Tensoes das 3 Fases
                             Vam2(k) ; Vbm2(k) ; Vcm2(k)];  %%Parte Imaginaria das Tensoes das 3 Fases
        end
    end
%%<

end               %%'end' para 'Vret' como condi��o inicial

%% Calculando os Valores de Refer�ncia para as Cargas ZIP
%%>
if Ent(5)==2 || Ent(5)==3
    Va0 = zeros(Ncz1,1);
    Vb0 = zeros(Ncz1,1);
    Vc0 = zeros(Ncz1,1);
    
    for i=1:Ncz
        Va0(i) = Va(find(kb==kcz(i)));
        Vb0(i) = Vb(find(kb==kcz(i)));
        Vc0(i) = Vc(find(kb==kcz(i)));
    end
    if Ent(1)==1 && Ncz2>0
        Va02 = zeros(Ncz2,1);
        Vb02 = zeros(Ncz2,1);
        Vc02 = zeros(Ncz2,1);
        for i=1:Ncz2
            Va02(i) = Va2(find(kb2==kcz2(i)));
            Vb02(i) = Vb2(find(kb2==kcz2(i)));
            Vc02(i) = Vc2(find(kb2==kcz2(i)));
        end
    end
end
%%<

% ======================================================================= %
% MONTAGEM DOS BLOCOS CONSTANTES (NAO ATUALIZADOS DURANTE O PROCESSO ITERATIVO)
% ======================================================================= %
%%OBS.: Esses blocos deveriam estar, de forma mais correta, dentro do Loop
%%Principal, onde seriam executados somente quando 'iter == 1'. Mas como �
%%praticamente imposs�vel entrar com as tens�es corretas e os c�lculos
%%abaixo s�o extremamente r�pidos, isso � feito externamente para melhorar
%%a visualiza��o do c�digo.

%%1>
%% MONTAGEM DA MATRIZ DE ADIMITANCIAS DO SISTEMA (Apenas para LT)
%%>
if Ent(3)==1
    %%A Primeira forma monta a matriz Completa, sem a Eliminacao das Barras
    %%de Referencia: Ela eh necessaria para o calculo dos Residuos
    Yret     = Ybarra_MICTH_7(h,Nb,Nl,Ne,p,q,00,lim,Raa,Rbb,Rcc,Rab,Rbc,Rca,Xaa,Xbb,Xcc,Xab,Xbc,Xca,Bash,Bbsh,Bcsh);   %%Permite a Inclus�o de 'Inf' para facilitar a remo��o de Linhas
    Yret_ref = Ybarra_MICTH_7(h,Nb,Nl,Ne,p,q,br,lim,Raa,Rbb,Rcc,Xaa,Xbb,Xcc,Bash,Bbsh,Bcsh);   %%Permite a Inclus�o de 'Inf' para facilitar a remo��o de Linhas
 	%Yret = Ybarra_MICTH_7(h,Nb,Nl,Ne,p,q,br,lim,Raa,Rbb,Rcc,Rab,Rbc,Rca,Xaa,Xbb,Xcc,Xab,Xbc,Xca,Bash,Bbsh,Bcsh);
	J0 = J0 + Yret_ref;
    Ybus=Yret;
end
%%<
% ======================================================================= %
% ======================================================================= %
%% INSER�AO DA CONTRIBUI�AO DOS TRANSFORMADORES EM TODAS AS FREQUENCIAS
%%>
if Ent(4)==1
    %%A Primeira forma monta a matriz Completa, sem a Eliminacao das Barras
    %%de Referencia: Ela eh necessaria para o calculo dos Residuos
	Yt     = Transf_4(h,Nb,Nt,Ne,pt,qt,nt,00,lim,Rtr,Xtr,Rf,Xm,Tap,yt);
    Yt_ref = Transf_4(h,Nb,Nt,Ne,pt,qt,nt,br,lim,Rtr,Xtr,Rf,Xm,Tap,yt);
    J0 = J0 + Yt_ref;
elseif Ent(4)==-1
    Yt     = Transf_4(h,Nb,Nt,Ne,pt,qt,[],00,lim,Rtr,Xtr,Rf,Xm,Tap,yt);
    Yt_ref = Transf_4(h,Nb,Nt,Ne,pt,qt,[],br,lim,Rtr,Xtr,Rf,Xm,Tap,yt);
    J0 = J0 + Yt_ref;
    Ybus=Ybus+Yt;
end
%%<
% ======================================================================= %
% ======================================================================= %
%% CONTRIBUICAO DAS CARGAS Zconst EM TODAS AS FREQUENCIAS
%%>
%%Para Todas as Freq. (Zconst: RLC)
if Ent(5)==1 || Ent(5)==3
    %%A Primeira forma monta a matriz Completa, sem a Eliminacao das Barras
    %%de Referencia: Ela eh necessaria para o calculo dos Residuos
	Yc     = JCargasRLC_h_4(w,h,Nb,Ne,kc,00,lim,RaC,RbC,RcC,CaC,CbC,CcC,LaC,LbC,LcC,y,c);
    Yc_ref = JCargasRLC_h_4(w,h,Nb,Ne,kc,br,lim,RaC,RbC,RcC,CaC,CbC,CcC,LaC,LbC,LcC,y,c);
    J0 = J0 + Yc_ref;
    Ybus=Ybus+Yc;
end
%%<
% ======================================================================= %
% ======================================================================= %
%% CONTRIBUICAO DAS CARGAS ZIP NAS FREQ. HARM�NICAS (CONFERIR!)
%%>
%%Para Todas as Freq. (Zconst: RLC)
if Ent(5)==2 || Ent(5)==3
    %%A Primeira forma monta a matriz Completa, sem a Eliminacao das Barras
    %%de Referencia: Ela eh necessaria para o calculo dos Residuos
	Ycz     = JCargas_h_6(F1,h,hf1,Nb,Ne,kcz,00,lim,Vret,Paz,Pbz,Pcz,Qaz,Qbz,Qcz,yz,zeros(Ncz1,1));
    Ycz_ref = JCargas_h_6(F1,h,hf1,Nb,Ne,kcz,br,lim,Vret,Paz,Pbz,Pcz,Qaz,Qbz,Qcz,yz,zeros(Ncz,1));
    
    if Ent(1)==1 && Ncz2>0
        Ycz2     = JCargas_h_6(F2,h,hf2,Nb,Ne,kcz2,00,lim,Vret,Paz2,Pbz2,Pcz2,Qaz2,Qbz2,Qcz2,yz2,zeros(Ncz2,1));
        Ycz2_ref = JCargas_h_6(F2,h,hf2,Nb,Ne,kcz2,br,lim,Vret,Paz2,Pbz2,Pcz2,Qaz2,Qbz2,Qcz2,yz2,zeros(Ncz2,1));
        
        Ycz      = Ycz     + Ycz2;
        Ycz_ref  = Ycz_ref + Ycz2_ref;
    end
    
    J0 = J0 + Ycz_ref;
    Ybus=Ybus+Ycz;
end
%%<
% ======================================================================= %
% ======================================================================= %
%% CONTRIBUICAO DOS BANCOS DE CAPACITORES DOS SVCs
%%>
if Ent(10)==1
    n   = length(Ras);   cbs = ones(n,1);                       %%R-C Serie
    Las = zeros(n,1);    Lbs = zeros(n,1);    Lcs = zeros(n,1);
    
	Yrs     = JCargasRLC_h_4(w,h,Nb,Ne,ks,00,lim,Ras,Rbs,Rcs,Cas,Cbs,Ccs,Las,Lbs,Lcs,ybs,cbs);
    Yrs_ref = JCargasRLC_h_4(w,h,Nb,Ne,ks,br,lim,Ras,Rbs,Rcs,Cas,Cbs,Ccs,Las,Lbs,Lcs,ybs,cbs);
    J0 = J0 + Yrs_ref;
    Ybus=Ybus+Yrs;
end, %Yrs = sparse(6*H*Nb + Ne , 6*H*Nb + Ne);
%%<
% ======================================================================= %
% ======================================================================= %
%% CONTRIBUICAO DOS GERADORES NAS FREQ. HARMONICAS
%%Logicamente nao existirah nenhum sistema sem um Gerador ao menos, mas para
%%tornar o codigo mais geral, irei fazer os testes para estes tambem.
%%>
if Ent(6)==1
    %%A Primeira forma monta a matriz Completa, sem a Eliminacao das Barras
    %%de Referencia: Ela eh necessaria para o calculo dos Residuos
	Yg     = JCargas_h_6(F1,h,hf1,Nb,Ne,kg,00,lim,Rg,Xg);
    Yg_ref = JCargas_h_6(F1,h,hf1,Nb,Ne,kg,br,lim,Rg,Xg);
	if Ent(1)==1
        Yg2     = JCargas_h_6(F2,h,hf2,Nb,Ne,kg2,00,lim,Rg2,Xg2);
        Yg2_ref = JCargas_h_6(F2,h,hf2,Nb,Ne,kg2,br,lim,Rg2,Xg2);
        Yg     = Yg     + Yg2;
        Yg_ref = Yg_ref + Yg2_ref;
	end
    J0 = J0 + Yg_ref;
    Ybus=Ybus+Yg;
end
%%<
% ======================================================================= %
% ======================================================================= %
%%1<

%%Passando Vret como Condi��o Inicial (Colocado aqui por causa das Cargas ZIP)
%>
if skip==0                              %%Sem Inicializa��o Pr�via
    %Vret = zeros(6*H*Nb,1);            %%Sem A��o nesse ponto do C�digo
elseif skip==0.5                        %%Com Inicializa��o Pr�via (Total ou Parcial)
    Vret = zeros(6*H*Nb,1);
    if length(Vret)>=length(V)
        Vret(1:length(V),1) = V;
    else
        Vret = V(1:length(Vret),1);
    end
end
%<

iter = 0;
end                                         %%'if' do skip inicial (Condicoes Iniciais)


%% Inicializa��o da Matriz J para os Conversores:
%%OBS.: Inicialmente isso seria mais eficiente, tanto do ponto de vista da
%%programacao, quanto do ponto de vista do esforco computacional total
%%necessario (existe um gasto de tempo para se criar essa matriz a cada
%%iteracao, caso esta seja formada por zeros). Entretanto, como passa a ser
%%necessario a existencia de duas dessas matrizes simultaneamente (uma aqui
%%e outra dentro de 'Conv_x'), para certas dimensoes isso se torna
%%impossivel, dado a memoria disponivel, o que forca a utilizacao direta de
%%matrizes esparsas. Isso por sua vez eleva o tempo necessario, uma vez que
%%utilizar uma matriz de zeros e posteriormente converte-la para 'sparse'
%%eh muito mais rapido. Por isso nao eh mais realizado!
%>
% if 6*H*Nb + Ne > Max
% %     Jj = sparse(6*H*Nb + Ne , 6*H*Nb + Ne);
%     Jj=0;                                 %%Apenas para compatibilizar
% else
% %     Jj = zeros(6*H*Nb + Ne , 6*H*Nb + Ne);
%     Jj=0;                                 %%Apenas para compatibilizar
% end
%<
% ======================================================================= %
%%                          LOOP PRINCIPAL
% ======================================================================= %
%%2>
while max(abs(DI))>Tol && iter <= cmax
    
    %% Relat�rio de N�o-Converg�ncia
    if iter>0
    %>
        aa = find(abs(DI(1:6*H*Nb))>Tol);
        a = [Barras(aa,1),Barras(aa,2)];
        a(a(:,1)==0,1)=Nb;
%         figure
% %           plot(a(:,2),a(:,1),'*','MarkerSize',7), grid on
%             plot(a(:,2),a(:,1),'s','MarkerFaceColor','k','MarkerSize',6), grid on
%             set(gca,'XTick',h(1):h(end))
%             title('Rela��o das Viola��es de Res�duos de Corrente')
%             xlabel('Harm�nicas')
%             ylabel('Barras')
%             zoom on
    end
    %<
    
    J = J0;
    
    if iter>0                               %%Para n�o calcular com as condi��es iniciais
    %%3>
    
% ======================================================================= %
% ======================================================================= %
%%  CONTRIBUICAO DOS GERADORES NA FREQ. FUNDAMENTAL
	%%>   Excluindo as Barras VTeta!
    if Ent(6)==1
		%>
        Jc = JCargas_6(F1,H,Nb,Ne,kg(xx1),br,Vret,PA(xx1),PB(xx1),PC(xx1),QA(xx1),QB(xx1),QC(xx1));
		J = J + Jc;
        
        if Ent(1)==1
            Jc = JCargas_6(F2,H,Nb,Ne,kg2(xx2),br,Vret,PA2(xx2),PB2(xx2),PC2(xx2),QA2(xx2),QB2(xx2),QC2(xx2));
			J = J + Jc;
        end
        %%fprintf('\n Ok: Geradores - ');
% 		clear Jc
    end
	%<
	%%<
% ======================================================================= %
% ======================================================================= %
%%  CONTRIBUICAO DAS CARGAS ZIP NA FREQ. FUNDAMENTAL (CONFERIR!)
	%%>   Excluindo as Barras VTeta!
    if Ent(5)==2 || Ent(5)==3
		%>
        Jc = JCargas_6(F1,H,Nb,Ne,kcz,br,Vret,Paz,Pbz,Pcz,Qaz,Qbz,Qcz,Pp,Pi,Pz,Qp,Qi,Qz,yz,Va0,Vb0,Vc0);
		J = J + Jc;
        
        if Ent(1)==1 && Ncz2>0
            Jc = JCargas_6(F2,H,Nb,Ne,kcz2,br,Vret,Paz2,Pbz2,Pcz2,Qaz2,Qbz2,Qcz2,pp2,pi2,pz2,qp2,qi2,qz2,yz2,Va02,Vb02,Vc02);
			J = J + Jc;
        end
        %%fprintf('Ok: Cargas ZIP - ');
% 		clear Jc
    end
	%<
	%%<
% ======================================================================= %
% ======================================================================= %
%%  INSER�AO DA CONTRIBUI�AO DOS RCTs (APENAS S�NCRONO)
	%%>
    if Ent(9)==1
        %>
        Jc = RCT_7(w,h,Nb,Ne,Nrct,krct,T0r_1,T0r_2,Tfr_1,Tfr_2,Rar,Rbr,Rcr,Lar,Lbr,Lcr,yrct,Max);
        %<
        J = J + Jc;
       % %fprintf('Ok: RCTs - ');
% 		clear Jc
    end
	%%<
% ======================================================================= %
% ======================================================================= %
%%  INSER�AO DA CONTRIBUI�AO DOS SVCs (APENAS S�NCRONO)
	%%>
    if Ent(10)==1
        %>
        Jc = RCT_7(w,h,Nb,Ne,Ns,ks,T0s_1,T0s_2,Tfs_1,Tfs_2,Rars,Rbrs,Rcrs,Lars,Lbrs,Lcrs,yrs,Max);
        %<
        J = J + Jc;
       % %fprintf('Ok: SCVs - ');
% 	    clear Jc
    end
    %%<
% ======================================================================= %
% ======================================================================= %
%%  INSER�AO DA CONTRIBUI�AO DOS RETIFICADORES
	%%>
    if Ent(7)==1
        %%Parte Superior:
        %>
%       Jc = Conv_6(w,h,hf2,Nb,Ne,Nr,kr,kr1,br,T0_1,Tf_1,Ra,Rb,Rc,L1,Max);
        Jc = Conv_7(w,h,hf2,Nb,Ne,Nr,kr,kr1,br,T0_1,Tf_1,Ra,Rb,Rc,L1,Max);
		%<
        J = J + Jc;
        %%Parte Inferior:
        %>
%       Jc = Conv_6(w,h,hf2,Nb,Ne,Nr,kr,kr2,br,T0_2,Tf_2,Ra,Rb,Rc,L2,Max);
        Jc = Conv_7(w,h,hf2,Nb,Ne,Nr,kr,kr2,br,T0_2,Tf_2,Ra,Rb,Rc,L2,Max);
        %<
        J = J + Jc;
        
        %fprintf('Ok: Retificadores - ');
% 		clear Jc
    end
    clc
	%%<
% ======================================================================= %
% ======================================================================= %
%%  INSER�AO DA CONTRIBUI�AO DOS INVERSORES
	%%>
    %%OBS.: Lembrar que as v�lvulas SUPERIORES Sempre disparam no
    %%Semi-Ciclo POSITIVO e as INFERIORES no NEGATIVO! (21/02/08)
    if Ent(8)==1
        %%Parte Superior:
        %>
%       Jc = Conv_6(w,h,hf1,Nb,Ne,Ni,ki,ki1,br,T0i_1,Tfi_1,Ria,Rib,Ric,Li1,Max);
        Jc = Conv_7(w,h,hf1,Nb,Ne,Ni,ki,ki1,br,T0i_1,Tfi_1,Ria,Rib,Ric,Li1,Max);
        %<
        J = J + Jc;
        %%Parte Inferior:
        %>
%       Jc = Conv_6(w,h,hf1,Nb,Ne,Ni,ki,ki2,br,T0i_2,Tfi_2,Ria,Rib,Ric,Li2,Max);
        Jc = Conv_7(w,h,hf1,Nb,Ne,Ni,ki,ki2,br,T0i_2,Tfi_2,Ria,Rib,Ric,Li2,Max);
        %<
        J = J + Jc;
        
        %fprintf('Ok: Inversores - ');
% 	    clear Jc
    end
    %%<
% ======================================================================= %
% ======================================================================= %
%%  EXPANSAO DO SISTEMA MATRICIAL (Eq. de Controle)
	%%>
	% % CONTROLE DE TENSAO POR REATIVOS (Apenas a Componente Fundamental):
    %%OBS.: As Tens�es de Sa�da da Fun��o abaixo n�o s�o necess�rias nesse ponto,
    %%%     mas por simplicidade s�o extra�das.
    if Ent(6)==1
		%>
		[VA,VB,VC,L_ad,C_ad] = JCont_PV_3(F1,H,Nb,kg,xx1,BarCtr,Vret,1);
        
        if Ent(1)==1
            [VA2,VB2,VC2,L_ad2,C_ad2] = JCont_PV_3(F2,H,Nb,kg2,xx2,BarCtr2,Vret,1);
            
            J(6*H*Nb+1:6*H*Nb+Ne , 1:6*H*Nb) = [L_ad ; L_ad2];
			J(1:6*H*Nb , 6*H*Nb+1:6*H*Nb+Ne) = [C_ad , C_ad2];
        else
            J(6*H*Nb+1:6*H*Nb+Ne , 1:6*H*Nb) = L_ad;
			J(1:6*H*Nb , 6*H*Nb+1:6*H*Nb+Ne) = C_ad;
        end
        
        %%fprintf('Ok: Eq. Controle - ');
    end
	%<
	%%<
% ======================================================================= %
% ======================================================================= %
%%  ELIMINA��O DAS BARRAS DE REFER�NCIA (VTeta): Frequ�ncia 1 (e 2)
    %%Essa parte foi feita internamente em todas as Rotinas, faltando
    %%apenas a colocacao de '1' nas diagonais
	%%>
    if Ent(6)==1
        for ir=1:length(br)
            for ih=1:H
                ii = 6*Nb*(ih - 1)  +  6*(br(ir)-1) + 1; %%Endere�o na Matriz 'J'
                
%               J(ii:ii+5,:) = zeros(6,6*H*Nb+Ne);  %%Zerando as Linhas
%               J(:,ii:ii+5) = zeros(6*H*Nb+Ne,6);  %%Zerando as Colunas
                
		        J(ii:ii+5,ii:ii+5) = eye(6);        %%'1' na Diagonal
                
% 	            DI(ii:ii+5,1) = zeros(6,1);         %%J� feito na Itera�ao anterior
            end
        end
        %fprintf('Ok: VTeta (VTeta_2) \n');
    end
    %%<
% ======================================================================= %
%%  ELIMINA��O DAS BARRAS DE REFER�NCIA (VTeta): Frequ�ncia 2
	%%>
    %%Incorporado no Bloco de Cima
    %%<
% ======================================================================= %
%%  ELIMINA��O DAS PARTES IMAGINARIAS DA FREQU�NCIA CC: TODAS AS BARRAS (NAO EXISTEM)
    %%Essa parte foi feita internamente em todas as Rotinas, faltando
    %%apenas a colocacao de '1' nas diagonais
	%%>
    for ih=1:H
        if h(ih)==0
            for ib=1:Nb
                ii = 6*Nb*(ih - 1)  +  6*(ib-1) + 1;%%Endere�o das posicoes Imaginarias
        
%               J(ii+3:ii+5,:) = zeros(3,6*H*Nb+Ne);%%Zerando as Linhas
%               J(:,ii+3:ii+5) = zeros(6*H*Nb+Ne,3);%%Zerando as Colunas
        
                J(ii+3:ii+5,ii+3:ii+5) = eye(3);    %%Colocando '1' na diagonal
        
%               DI(ii+3:ii+5,1) = zeros(3,1);       %%J� feito na Itera�ao anterior
            end
            break
        end
    end        
    %%<
% ======================================================================= %
%%  ELIMINA��O DAS POSI��ES INEXISTENTES PARA AS BARAS CC: TODAS AS BARRAS CC (Fases 'B' e 'C')
    %%Apenas a Coloca��o de '1' nas diagonais para as partes Reais. A
    %%matriz j� � montada com zeros para as outras posi��es e as partes
    %%Imagin�rias s�o feitas acima.
	%%>
    if iter==1
        k = kcc;                            %%Barras CC adicionais �s dos Conversores
        if Ent(7)==1
            k = [k ; kr1 ; kr2]; %#ok<AGROW>
        end
        if Ent(8)==1
            k = [k ; ki1 ; ki2]; %#ok<AGROW>
        end
%         k = k(:);                         %%Transformando em Vetor Coluna (N�o Precisa)
    end
    for ik=1:length(k)
        for ih=1:H
            i = 6*Nb*(ih - 1)  + 6*(k(ik)-1) + 1 + 1; %%Endere�o na Matriz 'J'

            if h(ih)==0
                J(i  :i+1,i  :i+1) = eye(2);        %%Acrescentando 1 (Apenas Real, pois o Imag. j� foi feito acima)
%               J(i  :i+4,i  :i+4) = eye(5);        %%Acrescentando 1 para as Fases 'B' e 'C' (Tudo)
            else
                J(i  :i+1,i  :i+1) = eye(2,2);      %%Acrescentando 1 na Diagonal (Imagin�rio)
                J(i+3:i+4,i+3:i+4) = eye(2,2);      %%Acrescentando 1 na Diagonal (Real)
            end
        end
    end
    %%<
% ======================================================================= %
%%  ELIMINA��O DAS LINHAS FICT�CIAS CRIADAS DO LADO CC: BARRAS DOS CONVERSORES (Feito acima)
%     %%Apenas a Coloca��o de '1' nas diagonais. A matriz j� � montada com
%     %%zeros para as outras posi��es.
% 	%%>
%     if iter==1
%         k = [];
%         if Ent(7)==1
%             k = [k , kr1, kr2]; %#ok<AGROW>
%         end
%         if Ent(8)==1
%             k = [k , ki1, ki2]; %#ok<AGROW>
%         end
%         k = k(:);                           %%Transformando em Vetor Coluna
%     end
%     if ~isempty(k)
%         for ik=1:length(k)
%             for ih=1:H
%                 i = 6*Nb*(ih - 1)  + 6*(k(ik)-1) + 2; %%Endere�o na Matriz 'J'
% 
%                 if h(ih)==0
% %                   J(i:i+4,:) = zeros(5,6*H*Nb+Ne);    %%Zerando as Linhas
% %                   J(:,i:i+4) = zeros(6*H*Nb+Ne,5);    %%Zerando as Colunas
% 
%                     J(i  :i+1,i  :i+1) = eye(2);        %%Acrescentando 1 (Apenas Real, pois o Imag. j� foi feito acima)
% %                   J(i  :i+4,i  :i+4) = eye(5);      %%Acrescentando 1 para CC (Tudo)
%                 else
% %                   J(i:i+1,:) = zeros(2,6*H*Nb+Ne);    %%Zerando as Linhas
% %                   J(:,i:i+1) = zeros(6*H*Nb+Ne,2);    %%Zerando as Colunas
% %                   J(i+3:i+4,:) = zeros(2,6*H*Nb+Ne);  %%Zerando as Linhas
% %                   J(:,i+3:i+4) = zeros(6*H*Nb+Ne,2);  %%Zerando as Colunas
% 
%                     J(i  :i+1,i  :i+1) = eye(2,2);      %%Acrescentando 1
%                     J(i+3:i+4,i+3:i+4) = eye(2,2);      %%Acrescentando 1
%                 end
%             end
%         end
%     end
%     %%<
% ======================================================================= %
%%  ELIMINA��O DAS POSI��ES INEXISTENTES PARA AS BARRAS DE NEUTRO (Fases 'A' ou 'B' e 'C')
    %%Apenas a Coloca��o de '1' nas diagonais. A matriz j� � montada com
    %%zeros para as outras posi��es.
	%%>
    if Ent(4)==1
        if iter==1
            [in,i,kn] = find(nt);           %%A segunda sa�da � Lixo!
        end
        Tr = zeros(1,length(kn));
        for ik=1:length(kn)
            In = in(ik);
            if ~sum(In==Tr)                 %%O ganho de tempo com isso � quase zero!
            
                if yt(In) == 00 || yt(In) == 10
                    sec = 0;
                else
                    %%Para aproveitar uma mesma Barra para Neutros de
                    %%Transformadores Diferentes (Um na Fase A e outro na Fase B)!!!
                    %>
                    in2 = find(nt==nt(In));
                    in2 = setdiff(in2,In);              %%Endere�o do outro Transf. que usa a Mesma Barra de Neutro
                    if ~isempty(in2)
%                   if yt(in2) == 00 || yt(in2) == 10
                    Tr(ik) = in2;
                        sec = 0;
                    else
                    %<
                        sec = 1;
                    end
                end
                if yt(In) == 00 || yt(In) == 01
                    pri = 0;
                else
                    %%Para aproveitar uma mesma Barra para Neutros de
                    %%Transformadores Diferentes!!!
                    %>
                    in2 = find(nt==nt(In));
                    in2 = setdiff(in2,In);              %%Endere�o do outro Transf. que usa a Mesma Barra de Neutro
                    if ~isempty(in2)
%                   if yt(in2) == 00 || yt(in2) == 01
                    Tr(ik) = in2;
                        pri = 0;
                    else
                    %<
                        pri = 1;
                    end
                end
                Diag = [sec 0 0 ; 0 pri 0 ; 0 0 1];     %%'0' indica nenhuma altera��o e '1' a Elimina��o de Linhas e Colunas
                for ih=1:H
                    i = 6*Nb*(ih - 1)  + 6*(kn(ik)-1) + 1;  %%Endere�o na Matriz 'J'

                    if h(ih)==0
                        %%Acrescentando 1 (Apenas Real, pois o Imag. j� foi feito acima)
                        J(i  :i+2,i  :i+2) = J(i  :i+2,i  :i+2) + Diag;
                    else
                        J(i  :i+2,i  :i+2) = J(i  :i+2,i  :i+2) + Diag; %%(Imagin�rio)
                        J(i+3:i+5,i+3:i+5) = J(i+3:i+5,i+3:i+5) + Diag; %%(Real)
                    end
                end         %%'for' do Jacobiano
            end             %%'if' da Barra de Neutro
        end                 %%'for' principal
    end
    %%<
% ======================================================================= %
% ======================================================================= %
%%  ATUALIZA��O DAS VARI�VEIS
    %%>
    %%OBS.: Utiliza os Res�duos de Corrente da Itera��o Anterior!!!

    DX = J\DI;                          %%Residuos de Tensao e Reativos

    Delta_V = DX(1:6*H*Nb,1);           %%Rediduos de Tensao
    Delta_Q = DX(6*H*Nb+1:end,1);       %%Residuos de Reativos

    Vret = Vret + Delta_V;              %%Atualizacao das Tensoes

    try                                 %#ok<TRYNC> %%Conferir!!!
        QA(xx1) = QA(xx1) + Delta_Q(1:3:3*Lxx1);     %#ok<AGROW> %%Atualizacao da Geracao de Reativos (Freq. 1)
        QB(xx1) = QB(xx1) + Delta_Q(2:3:3*Lxx1);     %#ok<AGROW> %%Atualizacao da Geracao de Reativos (Freq. 1)
        QC(xx1) = QC(xx1) + Delta_Q(3:3:3*Lxx1);     %#ok<AGROW> %%Atualizacao da Geracao de Reativos (Freq. 1)
        
        if Ent(1)==1
            QA2(xx2) = QA2(xx2) + Delta_Q(3*Lxx1+1:3:end);     %#ok<AGROW> %%Atualizacao da Geracao de Reativos (Freq. 2)
            QB2(xx2) = QB2(xx2) + Delta_Q(3*Lxx1+2:3:end);     %#ok<AGROW> %%Atualizacao da Geracao de Reativos (Freq. 2)
            QC2(xx2) = QC2(xx2) + Delta_Q(3*Lxx1+3:3:end);     %#ok<AGROW> %%Atualizacao da Geracao de Reativos (Freq. 2)
        end
    end
    
    %fprintf(' Ok: Atualiza��o \n\n');
    
    end         %%'if' para evitar a 1� itera��o (iter=0)
    %%3<
    
% ======================================================================= %
% ======================================================================= %
%%  ELIMINA��O DAS PARTES IMAGINARIAS DA COMPONENTE CC DA TENS�O (Apenas por Precau��o: Pode ser RETIRADO!!!)
%     %%N�o Faz sentido existir, mas como ainda n�o foi feita a limpeza
%     %%interna de Todas as Funcoes para a Freq CC (para alguns elementos na
%     %%Freq CC), essa elimina��o se faz neces�ria ==> NAO MAIS!!!.
% 	%%>
%     for ih=1:H
%         if h(ih)==0
%             for ib=1:Nb
%                 ii = 6*Nb*(ih - 1)  +  6*(ib-1) + 1;%%Endere�o das posicoes Imaginarias
%         
%                 Vret(ii+3:ii+5,1)  = zeros(3,1);
%             end
%             break
%         end
%     end        
%     %%<

% ======================================================================= %
    
% ====== x ===== x ===== x ===== x ===== x ===== x ===== x ===== x ====== %

% ======================================================================= %

%%  CALCULO DOS CRUZAMENTOS POR ZERO DOS RCTs E SVCs (TERMINAR!!!) (APENAS S�NCRONO)
	%%>
    if Ent(9)==1
        %%RCTs:
        %>
        [T0r_1,T0r_2,Tfr_1,Tfr_2] = CalcT_R_1(w,h,hv,hf1,Nb,Nrct,krct,yrct,Vret,Irct,...
                                                 Ton1,Toff1,Ton2,Toff2,iter,Iter0r);

        [Tfr_1 - T0r_1   (Tfr_1 - T0r_1)*120*180];
        [Tfr_2 - T0r_2   (Tfr_2 - T0r_2)*120*180];
        %<
    end
    if Ent(10)==1
        %%SVCs:
        %>
        [T0s_1,T0s_2,Tfs_1,Tfs_2] = CalcT_R_1(w,h,hv,hf1,Nb,Ns,ks,yrs,Vret,Isvc,...
                                              Ton1s,Toff1s,Ton2s,Toff2s,iter,Iter0r);
        %<
        [Tfs_1 - T0s_1   (Tfs_1 - T0s_1)*120*180];
        [Tfs_2 - T0s_2   (Tfs_2 - T0s_2)*120*180];
    end
    %%<
% ======================================================================= %
% ======================================================================= %
%%  CALCULO DOS CRUZAMENTOS POR ZERO DOS CONVERSORES (TERMINAR!!!) ==> CONFERIR!!!
	%%>
    %%OBS.: Esta originalmente associado 'hv2' e 'hf2' para o Retificador
    %%%     e 'hv' e 'hf1' para o Inversor. Assim, para trocar qual a
    %%%     frequencia configurada para cada, basta trocar no Arquivo de
    %%%     Dados (O Padrao eh 'hv' e 'hf1' para 60Hz e 'hv2' e 'hf2' para 50Hz)
    if Ent(7)==1
        %%Retificador Superior:
        %>
        [T0_1,Tf_1] = CalcT_2_1(w,h,hv2,hf2,Nb,Nr,kr,kr1,Vret,Ir1,1,1,alfa,Toff,...
                                                       Ra,Rb,Rc,L1,iter,Iter0);
        %<
        %%Retificador Inferior:
        %>
        [T0_2,Tf_2] = CalcT_2_1(w,h,hv2,hf2,Nb,Nr,kr,kr2,Vret,Ir2,2,1,alfa,Toff,...
                                                       Ra,Rb,Rc,L2,iter,Iter0);
        %<
        DTR1 = Tf_1 - T0_1;
        DTR2 = Tf_2 - T0_2;
    end
    if Ent(8)==1
        %%Inversor Superior:
        %%Lembrar que a refer�ncia eh sempre na Barra CA (CA - CC),
        %%independente de se tratar do Retificador ou do Inversor, por isso
        %%as Correntes tem sinais trocados.
        %>
		[T0i_1,Tfi_1] = CalcT_2_1(w,h,hv,hf1,Nb,Ni,ki,ki1,Vret,-Ii1,1,2,alfai,Toffi,...
                                                       Ria,Rib,Ric,Li1,iter,Iter0);
        %<
        %%Inversor Inferior:
        %>
        [T0i_2,Tfi_2] = CalcT_2_1(w,h,hv,hf1,Nb,Ni,ki,ki2,Vret,-Ii2,2,2,alfai,Toffi,...
                                                       Ria,Rib,Ric,Li2,iter,Iter0);
        %<
        DTI1 = Tfi_1 - T0i_1;
        DTI2 = Tfi_2 - T0i_2;
    end
    %%<
% ======================================================================= %
% ======================================================================= %
%% INICIALIZACAO DOS RESIDUOS DE CORRENTE
%>
DI  = zeros(6*H*Nb,1);                      %%As Eq de Controle Sao Adicionadas depois
%<
% ======================================================================= %
% ======================================================================= %
%%  CALCULO DOS RESIDUOS DE CORRENTE PARA OS GERADORES NA FREQUENCIA FUND.
	%%>   Excluindo as Barras VTeta!
    %>
    if Ent(6)==1
	% 	DIg1 = Res_PQ_3(F1,H,Nb,Ng,kg,Vret,PA,PB,PC,QA,QB,QC,1);    %%Modelo
        DIg1 = Res_PQ_4(F1,H,Nb,Lxx1,kg(xx1),Vret,PA(xx1),PB(xx1),PC(xx1),QA(xx1),QB(xx1),QC(xx1));
        
        if Ent(1)==1
            DIg1_2 = Res_PQ_4(F2,H,Nb,Lxx2,kg2(xx2),Vret,PA2(xx2),PB2(xx2),PC2(xx2),QA2(xx2),QB2(xx2),QC2(xx2));
            DIg1 = DIg1 + DIg1_2;
        end
        %<
% ======================================================================= %
% ======================================================================= %
%%  CALCULO DOS RESIDUOS DE CORRENTE PARA OS GERADORES NAS FREQ. HARMONICAS

		DIg2 = Res_Y(Ne,Yg,Vret);
        DIg = DIg1 + DIg2;                      %%Residuo Total
        
        DI = DI + DIg;
	end
	%%<
% ======================================================================= %
% ======================================================================= %
%%  CALCULO DOS RESIDUOS DE CORRENTE PARA AS CARGAS ZIP NA FREQUENCIA FUND. (CONFERIR!)
	%%>   Excluindo as Barras VTeta!
    %>
    if Ent(5)==2 || Ent(5)==3
        DIy1 = Res_PQ_4(F1,H,Nb,Ncz1,kcz,Vret,Paz,Pbz,Pcz,Qaz,Qbz,Qcz,Pp,Pi,Pz,Qp,Qi,Qz,Va0,Vb0,Vc0,yz);
        
        if Ent(1)==1 && Ncz2>0
            DIy1_2 = Res_PQ_4(F2,H,Nb,Ncz2,kcz2,Vret,Paz2,Pbz2,Pcz2,Qaz2,Qbz2,Qcz2,pp2,pi2,pz2,qp2,qi2,qz2,Va02,Vb02,Vc02,yz2);
            DIy1 = DIy1 + DIy1_2;
        end
        %<
% ======================================================================= %
% ======================================================================= %
%%  CALCULO DOS RESIDUOS DE CORRENTE PARA AS CARGAS ZIP NAS FREQ. HARMONICAS (CONFERIR!)

		DIy2 = Res_Y(Ne,Ycz,Vret);
        DIy = DIy1 + DIy2;                      %%Residuo Total
        
        DI = DI + DIy;
	end
	%%<
% ======================================================================= %
% ======================================================================= %
%%  CALCULO DOS RESIDUOS DE CORRENTE PARA A LINHA DE TRANSMISSAO
    %%>
    if Ent(3)==1
        DIl = Res_Y(Ne,Yret,Vret);
        
        DI = DI + DIl;
    end
    %%<
% ======================================================================= %
% ======================================================================= %
%%  CALCULO DOS RESIDUOS DE CORRENTE PARA OS TRANSFORMADORES
    %%>
    if Ent(4)==1 || Ent(4)==-1
        DIt = Res_Y(Ne,Yt,Vret);
        
        DI = DI + DIt;
    end
    %%<
% ======================================================================= %
% ======================================================================= %
%%  CALCULO DOS RESIDUOS DE CORRENTE PARA AS CARGAS (RLC)
	%%>
    if Ent(5)==1 || Ent(5)==3
        DIy0 = Res_Y(Ne,Yc,Vret);
        if Ent(5)==3
            DIy = DIy + DIy0;
        else
            DIy = DIy0;
        end        
        DI = DI + DIy0;
    end
	%%<
% ======================================================================= %
% ======================================================================= %
%%  CALCULO DOS RESIDUOS DE CORRENTE PARA AS FONTES DE CORRENTES HARM�NICAS
	%%>
    if Ent(12)==1
        if iter==0
            for i=1:Ninj
                DInj(6*Nb*(hinj(i)-1) +1 + kinj(i)*6-6) = Inja(i,1);  %%Real A
                DInj(6*Nb*(hinj(i)-1) +1 + kinj(i)*6-3) = Inja(i,2);  %%Imag A
                
                DInj(6*Nb*(hinj(i)-1) +1 + kinj(i)*6-5) = Injb(i,1);  %%Real B
                DInj(6*Nb*(hinj(i)-1) +1 + kinj(i)*6-2) = Injb(i,2);  %%Imag B
                
                DInj(6*Nb*(hinj(i)-1) +1 + kinj(i)*6-4) = Injc(i,1);  %%Real C
                DInj(6*Nb*(hinj(i)-1) +1 + kinj(i)*6-1) = Injc(i,2);  %%Imag C
            end
        end

        DI = DI + DInj;
    end
	%%<
% ======================================================================= %
% ======================================================================= %
%%  CALCULO DOS RESIDUOS DE CORRENTE PARA OS RCTs  (APENAS S�NCRONO)
	%%>
    if Ent(9)==1
        %>
        [DIrct,Irct] = Res_RCT_5(w,h,Nb,Nrct,krct,Vret,T0r_1,T0r_2,Tfr_1,Tfr_2,Rar,Rbr,Rcr,Lar,Lbr,Lcr,yrct);
        %<
        
        DI = DI + DIrct;
        %fprintf('\n Ok: Res. RCTs - ');
    end
	%%<
% ======================================================================= %
% ======================================================================= %
%%  CALCULO DOS RESIDUOS DE CORRENTE PARA OS SVCs (APENAS S�NCRONO) (CONFERIR!!!)
	%%>
    if Ent(10)==1
        %>
        [DIsvc,Isvc] = Res_RCT_5(w,h,Nb,Ns,ks,Vret,T0s_1,T0s_2,Tfs_1,Tfs_2,Rars,Rbrs,Rcrs,Lars,Lbrs,Lcrs,yrs);
        DIbs  = Res_Y(Ne,Yrs,Vret);             %%Para o Banco de Capacitores
        DIsvc = DIsvc + DIbs;
         %<
        
        DI = DI + DIsvc;
        %fprintf('Ok: Res. SVCs - ');
    end
    %%<
% ======================================================================= %
% ======================================================================= %
%%  CALCULO DOS RESIDUOS DE CORRENTE PARA OS RETIFICADORES
	%%>
    if Ent(7)==1
        %%Parte Superior:
        %>
        [DI1,Ir1] = Res_Conv_5(w,h,hf2,Nb,Nr,kr,kr1,Vret,T0_1,Tf_1,Ra,Rb,Rc,L1,1,1);
        %<
        %%Parte Inferior:
        %>
        [DI2,Ir2] = Res_Conv_5(w,h,hf2,Nb,Nr,kr,kr2,Vret,T0_2,Tf_2,Ra,Rb,Rc,L2,2,1);
        %<
        DIr = DI1 + DI2;
        %<
        
        if Ent(1)==1         %%Integrando, dentro de 1 periodo (F.S.), todos os interalos de conducao
            T01 = T0_1;        Tf1 = Tf_1;
            T02 = T0_2;        Tf2 = Tf_2;
            
            for i=2:hf2
                T01 = T01 + 2*pi/(w*hf2);            Tf1 = Tf1 + 2*pi/(w*hf2);
                T02 = T02 + 2*pi/(w*hf2);            Tf2 = Tf2 + 2*pi/(w*hf2);
                
                [DI1,Ir] = Res_Conv_5(w,h,hf2,Nb,Nr,kr,kr1,Vret,T01,Tf1,Ra,Rb,Rc,L1,1,1);
                DIr = DIr + DI1;
                Ir1 = Ir1 + Ir;
                [DI2,Ir] = Res_Conv_5(w,h,hf2,Nb,Nr,kr,kr2,Vret,T02,Tf2,Ra,Rb,Rc,L2,2,1);
                DIr = DIr + DI2;
                Ir2 = Ir2 + Ir;
            end
        end
        DI = DI + DIr;
        %fprintf('Ok: Res. Retificadores - ');
    end
	%%<
% ======================================================================= %
% ======================================================================= %
%%  CALCULO DOS RESIDUOS DE CORRENTE PARA OS INVERSORES
    %%>
    if Ent(8)==1
        %%Parte Superior:
        %>
        [DI1,Ii1] = Res_Conv_5(w,h,hf1,Nb,Ni,ki,ki1,Vret,T0i_1,Tfi_1,Ria,Rib,Ric,Li1,1,2);
        %<
        %%Parte Inferior:
        %>
        [DI2,Ii2] = Res_Conv_5(w,h,hf1,Nb,Ni,ki,ki2,Vret,T0i_2,Tfi_2,Ria,Rib,Ric,Li2,2,2);
        %<
        DIi = DI1 + DI2;
        %<
        
        if Ent(1)==1         %%Integrando, dentro de 1 periodo (F.S.), todos os interalos de conducao
            T01 = T0i_1;        Tf1 = Tfi_1;
            T02 = T0i_2;        Tf2 = Tfi_2;
            
            for i=2:hf1
                T01 = T01 + 2*pi/(w*hf1);            Tf1 = Tf1 + 2*pi/(w*hf1);
                T02 = T02 + 2*pi/(w*hf1);            Tf2 = Tf2 + 2*pi/(w*hf1);
                
                [DI1,Ii] = Res_Conv_5(w,h,hf1,Nb,Ni,ki,ki1,Vret,T01,Tf1,Ria,Rib,Ric,Li1,1,2);
                DIi = DIi + DI1;
                Ii1 = Ii1 + Ii;
                [DI2,Ii] = Res_Conv_5(w,h,hf1,Nb,Ni,ki,ki2,Vret,T02,Tf2,Ria,Rib,Ric,Li2,2,2);
                DIi = DIi + DI2;
                Ii2 = Ii2 + Ii;
            end
        end
        DI = DI + DIi;
        %fprintf('Ok: Res. Inversores \n');
    end
	%%<
% ======================================================================= %
% ======================================================================= %
%%  TESN�ES NAS BARRAS CONTROLADAS
	%%>
    if Ent(6)==1
		% % CONTROLE DE TENSAO POR REATIVOS (Apenas a Componente Fundamental):
		%>
        [VA,VB,VC] = JCont_PV_3(F1,H,Nb,kg,xx1,BarCtr,Vret,0);
        
        if Ent(1)==1
            [VA2,VB2,VC2] = JCont_PV_3(F2,H,Nb,kg2,xx2,BarCtr2,Vret,0);
        end
		%<
% ======================================================================= %
% ======================================================================= %
%%  CALCULO DOS RESIDUOS DE TENS�O (Eq. de Controle)
	%>
        for i=1:Lxx1                            %%CONFERIR aqui!
            ii = 3*i-2;
            
            DV(ii:ii+2) = [VCtr(xx1(i))-VA(i) ; VCtr(xx1(i))-VB(i) ; VCtr(xx1(i))-VC(i)];
        end
        
        ix = 0;
        if Ent(1)==1
            for i=Lxx1+1:Lxx                    %%CONFERIR aqui!
                ix = ix+1;
                ii = 3*i-2;
                
                DV(ii:ii+2) = [VCtr2(xx2(ix))-VA2(ix) ; VCtr2(xx2(ix))-VB2(ix) ; VCtr2(xx2(ix))-VC2(ix)];
            end
        end
        %<
        DI = [DI ; DV];                     %#ok<AGROW>%%� Insignificante o tempo gasto
	end
	%%<
% ======================================================================= %
% ======================================================================= %
%%  SOMA DOS RES�DUOS (Nao Existe mais. Eh Feita por Partes)
%     %%>
% 	DI = [(DIg + DIl + DIy + DIr + DIi + DIt) ; DV];    %%Residuos de Corrente (Completo)
%     %%<
% ======================================================================= %
% ======================================================================= %
%%  ELIMINA��O DOS RES�DUOS DAS BARRAS DE REFER�NCIA (VTeta): Frequ�ncia 1 (e 2)
	%%>
    if Ent(6)==1
        for ir=1:length(br)
            for ih=1:H
                ii = 6*Nb*(ih - 1)  +  6*(br(ir)-1) + 1; %%Endere�o na Matriz 'J'
                
                DI(ii:ii+5,1) = zeros(6,1);
            end
        end
% ======================================================================= %
%%  ELIMINA��O DOS RES�DUOS DAS BARRAS DE REFER�NCIA (VTeta): Frequ�ncia 2
    %%>
    %%Incorporado no Bloco de Cima
    %%<
    end
    %%<
% ======================================================================= %
%%  ELIMINA��O DAS PARTES IMAGINARIAS DOS RES�DUOS DA COMPONENTE CC (Pode ser RETIRADO!!!)
%     %%Acho que essa parte podera ser posteriormente retirada (Quando a limpeza
%     %%do codigo for realizada)
% 	%%>
%     for ih=1:H
%         if h(ih)==0
%             for ib=1:Nb
%                 ii = 6*Nb*(ih - 1)  +  6*(ib-1) + 1;%%Endere�o das posicoes Imaginarias
%         
%                 DI(ii+3:ii+5,1)  = zeros(3,1);
%                 DIt(ii+3:ii+5,1) = zeros(3,1);      %%Para os Transf tb
%             end
%             break
%         end
%     end
%     %%<
% ======================================================================= %
%%  ELIMINA��O DAS POSI��ES INEXISTENTES PARA AS BARAS CC (Fases 'B' e 'C')
%%Acho que essa parte podera ser posteriormente retirada (Quando a limpeza
%%do codigo for realizada)
	%%>
    
    %%<
% ======================================================================= %
%%  ELIMINA��O DAS POSI��ES INEXISTENTES PARA AS BARRAS DE NEUTRO (Fases 'A' ou 'B' e 'C')
%%Acho que essa parte podera ser posteriormente retirada (Quando a limpeza
%%do codigo for realizada)
	%%>
    
    %%<
    
%% ===================================================================== %%
% ======================================================================= %
% ======================================================================= %
    
    Vdata(:,iter+1) = Vret;                 %%Historico de Tensoes
    
%   iter                                    %#ok<NOPTS> %%Na primeira vez que passa por aqui ainda n�o realizou nenhuma itera��o! (S�o os dados informados).
    iter = iter + 1;
    
    %beep
    Ti = tempo;
    tempo = toc;                            %%Tempo Acumulado
    
%     disp('  ');
%     disp('-----------------------------------------------------------');
%     %fprintf(' Itera��es Realizadas:\t\t\t\t\t\t%g',iter-1);
%     %fprintf('\n Tempo Gasto Itera��o/Simula��o:\t\t\t%g/%g',tempo-Ti,tempo+tempo1);
%     %fprintf('\n M�ximo Res�duo encontrado at� o momento:\t%g\n',max(abs(DI)));
%     disp('-----------------------------------------------------------');
    
    if exist('Nome','var')
        Jc = 0;
        save (Nome)                         %%Se Informado Salva o Workspace a cada itera��o
    end
end
%%2<    Fim do LOOP Principal

if skip==1
    tempo = tempo + tempo1;
end

iter = iter-1;                              %%Apenas para compatibilizar

Vdata = full(Vdata(:,1:iter+1));            %%Eliminando as colunas nao utilizadas


% ======================================================================= %
% =====================   xxx  xxx  xxx  xxx  xxx   ===================== %
% ======================================================================= %
%% FIM!!!
toc
Sp = nnz(J)/numel(J)*100;
% %fprintf('Esparsidade do Sistema: %g%%\n',100-Sp);
% Hora = clock;
% %fprintf('\n Simula��o Terminada:');
% %fprintf('\n Dia %02.0f/%02.0f/%4.0f �s %02.0f:%02.0f:%05.2f \n \n',Hora(3),Hora(2),Hora(1),Hora(4),Hora(5),Hora(6));
% load Som_Fim;                               %%Original em 22050 Hz
% wavplay(Som_Fim,8e3)

% clear Jc J0 Som_Fim

% % Apenas para facilitar:
% k_elo = unique([kr;kr1;kr2;ki;ki1;ki2]);                        %%Barras do Elo (6 ou 12 Pulsos)
% k_trf = unique([full(kT(2));pt(pt>full(kT(14)))]);              %%Barras do Elo com Transf. (6 ou 12 Pulsos)
% % ITotal = Ir1+Ir2+Ii1+Ii2;                                       %%Apenas Correntes no ELO (6 Pulsos)
% Id = 'pt';


% if Ent(1)==1
%     F = F2; %%Apenas para Compatibilizar
% 	Plot_Data_1
% end
% 
% F = F1; %%Apenas para Compatibilizar
% Plot_Data_1


% Saida_Dados_1

% profile off
% profile viewer
% % profile report