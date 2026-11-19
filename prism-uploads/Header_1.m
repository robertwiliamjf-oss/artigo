
%-------------------------------------------------------------------%
% % %        Entrada de Dados        % % %
%-------------------------------------------------------------------%

% clc
% lambda=1;

%% Inicializa��es:
skip = 0;                                   %%'1'   Entra com todas as condi��es iniciais (pode ser utilizado para retomar uma simulacao interrompida),
                                            %%'0'   Roda Normalmente
if skip==0
    clc
    skip=0;
elseif skip==0.5                            %%Inicializa��o Total ou Parcial
    clc
    skip;
    try
        V = Vret;
%         clearvars -except V skip
    catch                                   %#ok<CTCH>
        skip=0;
    end
end

if skip~=1
%%>
Tol = 1e-5;                                 %%Tolerancia (Padr�o)
cmax = 12;                                  %%Numero maximo de itera�oes permitidas
Max = 5e3;                                  %%M�xima dimens�o para se usar 'zeros' a 'sparse'
lim = 1e10;                                 %%Valor 'infinito' a ser adicionado:
Iter0 = cmax+10;                            %%Itera��o a partir da qual se calculam os cruzamentos por zero da Corrente nos Conversores
Iter0r= cmax+10;                            %%Itera��o a partir da qual se calculam os cruzamentos por zero da Corrente nos RCTs e SVCs
%%<

%% Carregando os Dados
%%>
% ======================================================================= %
x=1;                                        %%Tipo de Arquivo de dados a ser utilizado ('1' ou '5')
if x==1                                     %%Sistema S�NCRONO
%% S�ncrono

    
     Dados33_svc;             %%14barras

     %Dados; % 33 barras

elseif x==5                                 %%Interconexao ASSINCRONA

end


%% Compatibiliza��o
%>
ABC = 1;                                    %%ABC: '1' (Sequ�ncia ABC escrita Correta); ABC(troca): '0'; ABC(-): '-1' (Arquivos ANTIGOS COM YD)
%%SEPOPE: ABC = 0;

try
    Ent;
catch %#ok<CTCH>                            %%Para Bancos de Dados Antigos
    if x==1, Ass=0; elseif x==5, Ass=1; else error('Valor Errado para "x"! Deve ser "1" ou "5"!'), end
%           1		2		3		4		5		6		7		8		9		10		11		12.
%           Ass		Bar		Lin		Transf	Carg	Ger		Ret		Inv           RCT		SVC		Bcc		iinj
    Ent = [	Ass		1		1		-1		1		1		1		1		0		0		1		1];
end
%<
end


%% Fun��o Principal:

MICTH;
%% Sa�da (solu��o)
DI_Calc = -[DIl];
hh=1; % vai imprimir resultados para a ordem harmonica hh
[V,T,P,Q,I,Ti] = Saida_Dados_h_4(w,kB,h,hh,Nb,Tol,iter,cmax,Vret,DI,DI_Calc,0,0,0,1,0);
%[V,T,P,Q,I,Ti] = Saida_Dados_h_4(w,kB,h,hh,Nb,Tol,iter,cmax,Vret,DI,DI_Calc,0,0,0,1,0);

% Para obter as tens�es s� para uma frequ�ncia:
T = T*pi/180;  % Transforma��o de angulo em graus
Ti = Ti*pi/180;  % Transforma��o de angulo em graus




%% Resultados fasoriais por ordem harm�nica (fator)
fator=1; % qual ordem harm�nica deseja ver fasores

if fator~=1
    for i=1:(fator-1)
V([1:Nb],:)=[];
T([1:Nb],:)=[];
I([1:Nb],:)=[];
Ti([1:Nb],:)=[];
    end
end


% Tens�es fasoiais
for i=1:Nb
    [a,b]=pol2cart(T(i,3),V(i,3));
    Vaa(1,i) = a + 1i*b;

    [a,b]=pol2cart(T(i,4),V(i,4));
    Vbb(1,i) = a + 1i*b;

    [a,b]=pol2cart(T(i,5),V(i,5));
    Vcc(1,i) = a + 1i*b;
  
end

% Correntes injetadas fasoriais
for i=1:Nb
    [a,b]=pol2cart(Ti(i,3),I(i,3));
    Iaa(1,i) = a + 1i*b;

    [a,b]=pol2cart(Ti(i,4),I(i,4));
    Ibb(1,i) = a + 1i*b;

    [a,b]=pol2cart(Ti(i,5),I(i,5));
    Icc(1,i) = a + 1i*b;
  
end

 
 









% [THD,THD2] = THD_5(kB,1,3,[],[],h,hh,Nb,0,V,1,'pt',1)





% Yret(:,[19 20 21]) = [];

% THD:
% 
Id='pt';

global color                                            %#ok<REDEF>
color = [0.3 0.3 0.3 ; 0.6 0.6 0.6 ; 0.9 0.9 0.9];      %#ok<NASGU> %%Mapa de cores em Escala de Cinza
kL=[];
%%Identificacao da Freq. Fundamental
%>
kf = zeros(Nb,1);
kf(kb) = hf1;                               %%Barras com 'Freq 1'
try                                         %#ok<TRYNC>
	kf(kb2)= hf2;                           %%Barras com 'Freq 2'
end

Kk = full(kT(1:Nb));                        %%Barras para THD (Tensao)
hf = kf(Kk);                                %%N�mero da Harm�nica correspondente � Fundamental para as Barras desejadas (Tensao)
hfl= kf(kL);                                %%N�mero da Harm�nica correspondente � Fundamental para as Barras com Tens�o de Linha
%<

[THDv,THD2] = THD_5(kB,Kk,hf,[],hfl,h,1:length(h),Nb,1e-3,Vret,1,Id,0);                    %%Tens�es (Fase e Linha para 12 Pulsos)

Kk = full(kT(1:Nb));                       %%Barras para THD (Corrente)
hf = kf(Kk);                                %%N�mero da Harm�nica correspondente � Fundamental para as Barras desejadas (Corrente)
THDi        = THD_5(kB,Kk,hf,[],[],h,1:length(h),Nb,1e-3,-DIl,2,Id,0);                    %%Correntes no HVDC



