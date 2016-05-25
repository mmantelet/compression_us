function data=traitement_compression(rep,ws)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Description                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fonction permettant de r�cup�rer les donn�es US (Temps de vol, SOS,
% amplitudes, fr�quences,...) et de les synchroniser avec les donn�es du
% TAXT (position, force, temps de compression).
% % % % % % % % % % % % % % % % Entr�es % % % % % % % % % % % % % % % % % %
% rep = nom du repertoire o� se trouve le fichier data.mat comportant A :
% les donn�es brutes US et B : les donn�es rutes du TAXT.
% ws = largeur de la fen�tre utilis� pour le moyennage glissant.
% % % % % % % % % % % % % % % % Sorties % % % % % % % % % % % % % % % % % %
% data est une structure comportant toutes les donn�es trait�es et
% synchronis�es :
% tz(1ou 2) : position temporelle des zeros crossing de l'echo (1 ou 2)
% (nSignal X nzeros crossing)
% tm(1 ou 2) : position temporelle des maximums de l'echo (1 ou 2)
% (nSignal X n maximums)
% Az(1 ou 2) : amplitudes des zeros crossing de l'echo (1 ou 2)
% (nSignal X nzeros crossing)
% Am(1 ou 2) : amplitudes des maximums crossing de l'echo (1 ou 2)
% (nSignal X n maximums)
% P : Position synchronis�e en mm
% (n Signal X 1)
% F : Force Synchronis�e en N
% (n Signal X 1)
%  Sigf : signaux US filtr�es
% (n Signal X taille Signal)
%  Sig : signaux US
% (n Signal X taille Signal)
%  SigH : transform�e de Hilbert des signaux
% (n Signal X taille Signal)
%  SigHf : transform�e de Hilbert filtr�es des signaux
% (n Signal X taille Signal)
% spectreS(1 ou 2) : Spectre des echos (1 ou 2)
% (n Signal X taille spectre)
% spectrefpic(1 ou 2) : fr�quence du maximums des spectres des echos (1 ou 2)
% (n Signal X 1)
% spectrefpic(1 ou 2) : amplitude du maximum des spectres des echos (1 ou 2)
% (n Signal X 2)
% TOF : Temps de vol pour chaque signaux en secondes
% (n Signal X n zeros crossing )
% SOS : vitesse calcul�es � partir des TOF en m.s-1
% (n Signal X n zeros crossing )
% rAMP : Att�nuations calcul�es entre les maximums de l'�cho1 et ceux de
% l'�cho 2
% (n Signal X n Maximums)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rep='comp_5_020516';
% ws=5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Importation / Calibrage fen�tre d'int�r�t / Traitement du bruit      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([rep,'\data'])


Sig=(A(1000:end,:))';
[nS,lS]=size(Sig);      % nS: Nombre de signaux // lS: longueur des signaux
T_us=A(2,:)';           % T_us: Base de temps des �missions de signaux US          
T_us=T_us-T_us(1);

% Filtrage du bruit avec moyenne glissante
b=1/ws*ones(1,ws);      % Pond�ration du moyennage
a=1;                    % Valeur par d�faut dans la fonction
for i=1:nS
    y=filter(b,a,Sig(i,:));
    Sigf(i,:)=Sig(i,:);
    Sigf(i,1:end-round((ws-1)/2))=y(1+round((ws-1)/2):end);
    SigH(i,:)=abs(hilbert(Sig(i,:)));
    SigHf(i,:)=abs(hilbert(Sigf(i,:)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      Traitement du 1er signal                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculs de Temps de vol pour le signal de r�f�rence
F=1E6;                  %Fr�quence du pulse
T=1/F*1*(100E6);        %P�riode temporelle du signal
S=Sigf(1,:);            %Signal Filtr�
H=SigHf(1,:);           %Transform�e de Hilbert du signal filtr�
nP=5;                   %Nombre de pics � d�tecter
nE=2;                   %Nombre d'�chos analys�s
seuil1=0.4;             %hauteur minimal de d�tection des pics sur la transform�e de Hilbert. Param�tre compris entre 0 et 1
seuil2=200;             %�cart minimal entre 2 pics sur la transform�e de Hilbert
coeff_nbp=1;            %coefficient multiplicateur du nombre de point du signal (en cas d'interpolation coeff_nbp >1)
echo_fen=[100, 200];    %coordonn�e de la fen�tre d'un echo par rapport � zc1 et zc4. ex de coord [zc1-echo_fen(1), zc4+echo_fen(2)]

[echo_win, tm, tz, ~, ~, pksH,locsH,widths]=extract_tcarac_ver10(S, seuil1, seuil2, 3*T, nP, coeff_nbp, 0);


tfE(1)=echo_win(1,2)+40;    %fin de l'echo1
tfE(2)=echo_win(2,1)-40;    %d�but de l'�cho2

% D�finition du signal de r�f�rence ( fin de l'�cho 1 jusqu'au d�but de
% l'�cho 2)
Sigref=zeros(1,lS);
Sigref(1,tfE(1):tfE(2))=mean(Sigf(1:50,tfE(1):tfE(2)));

% Soustraction du signal de r�f�rence au reste de l'ensemble des signaux
for i=1:nS
    Sigf(i,:)=Sigf(i,:)-Sigref(1,:);
    SigHf(i,:)=abs(hilbert(Sigf(i,:)));
end


tm1(1,:)=tm.echo(1).max(:,1)';
tm2(1,:)=tm.echo(2).max(:,1)';
Am1(1,:)=tm.echo(1).max(:,2)';
Am2(1,:)=tm.echo(2).max(:,2)';

 
tz1(1,:)=tz.echo(1).zc(:,1)';
tz2(1,:)=tz.echo(2).zc(:,1)';
Az1(1,:)=tz.echo(1).zc(:,2)';
Az2(1,:)=tz.echo(2).zc(:,2)';

% tm : num de l'echo X num du temps du max
% Am : num de l'�cho X num de l'amplitude du max


[fpic1(1,:), A6dB1(1,:), apic1(1,:), Dfreq1]=freq_pulse_ver2(S, [tz1(1,1)-echo_fen(1) tz1(1,4)+echo_fen(2)], coeff_nbp, 1);
spectreS1(1,:)=Dfreq1.echo.fftS;
f=Dfreq1.echo.f;

[fpic2(1,:), A6dB2(1,:), apic2(1,:), Dfreq2]=freq_pulse_ver2(S, [tz2(1,1)-echo_fen(1) tz2(1,4)+echo_fen(2)], coeff_nbp, 1);
spectreS2(1,:)=Dfreq2.echo.fftS;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Traitement de tous les autres signaux de la s�lection            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=2:nS
    
    [tm1n, tz1n]=suivi_sig(Sigf(i,:), tm1(i-1,:), tz1(i-1,:), 0.2*T);
    
    tm1(i,:)=tm1n(:,1)';
    Am1(i,:)=tm1n(:,2)';
    
    tz1(i,:)=tz1n(:,1)';
    Az1(i,:)=tz1n(:,2)';
    
    [fpic1(i,:), A6dB1(i,:), apic1(i,:), Dfreq1]=freq_pulse_ver2(Sigf(i,:), [tz1(i,1)-echo_fen(1) tz1(i,4)+echo_fen(2)], coeff_nbp, 0);
    spectreS1(i,:)=Dfreq1.echo.fftS;

    
    [tm2n, tz2n]=suivi_sig(Sigf(i,:), tm2(i-1,:), tz2(i-1,:), 0.2*T);
    
    tm2(i,:)=tm2n(:,1)';
    Am2(i,:)=tm2n(:,2)';
    
    tz2(i,:)=tz2n(:,1)';
    Az2(i,:)=tz2n(:,2)';
    
    [fpic2(i,:), A6dB2(i,:), apic2(i,:), Dfreq2]=freq_pulse_ver2(Sigf(i,:), [tz2(i,1)-echo_fen(1) tz2(i,4)+echo_fen(2)], coeff_nbp, 0);
    spectreS2(i,:)=Dfreq2.echo.fftS;

    
    clear tm1n tz1n tm2n tz2n Dfreq1 Dfreq2

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               Chargement des donn�es du TAXTPlus                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T_co=B(:,3);          % Temps
P=B(:,2);             % Position
P=10-P;
F=B(:,1);             % Force


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Synchronisation et d�finitions des donn�es finales               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% d�tectection du point de synchro TAXT

[n1, ~]=detect_comp(P);
T_co=T_co-T_co(n1);
T_co=1e-3*round(1e3*T_co);

%% D�tection du point de synchro TOF

[debut_us, ~]=detect_comp(tz2(:,1));
T_us=T_us-T_us(debut_us);%D�finition de l'origine des temps US
T_us=1e-3*round(1e3*T_us);

%% Synchronisation entre les deux types de donn�es

t_1=max(min(T_us),min(T_co));   %temps de la premi�re acquisition du dernier acquisiteur � etre lanc�
t_2=min(max(T_us),max(T_co));   %temps de la derni�re acquisition du dernier acquisiteur � etre �teint

I=find(T_us>=t_1&T_us<=t_2);    %recherche du vecteur du vecteur interpol� sur les temps US
tf=T_us(I);                     %vecteur des temps US variant dans la m�me gamme que le vecteur interpol�

Pf=interp1(T_co,P,tf)';         %interpolation de la position en fonction du temps avec le vecteur des temps US
Ff=interp1(T_co,F,tf)';         %interpolation de la force en fonction du temps avec le vecteur des temps US

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%                  Mise en forme des donn�es primaires                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

data.temps=tf;
data.tz1=tz1(I,:);
data.tz2=tz2(I,:);
data.Az1=Az1(I,:);
data.Az2=Az2(I,:);
data.tm1=tm1(I,:);
data.tm2=tm2(I,:);
data.Am1=Am1(I,:);
data.Am2=Am2(I,:);
data.P=Pf;
data.F=Ff;
data.Sigf=Sigf(I,:);
data.Sig=Sig(I,:);
data.SigH=SigH(I,:);
data.SigHf=SigHf(I,:);
data.spectreS1=spectreS1;
data.spectreS2=spectreS2;
data.spectrefpic1=fpic1(I,:);
data.spectrefpic2=fpic2(I,:);
data.spectreapic1=apic1(I,:);
data.spectreapic2=apic2(I,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                          Calcul vitesse                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

[nS,lS]=size(data.Sigf);
SOS=zeros(nS,nP);
for i=1:nS
    for j=1:length(tz2(1,:))
        TOF(i,j)=((data.tz2(i,j)-data.tz1(i,j))*1e-8);
        SOS(i,j)=2*data.P(i)*1e-3./TOF(i,j);
        SOSf(:,j)=filter(b,a,SOS(:,j));
        SOS(1+ws:end,j)=SOSf(1+ws:end,j);
    end
    for j=1:length(tm2(1,:))
        rAMP(i,j)=20*log10(abs(data.Am2(i,j)/data.Am1(i,j)))/(2*0.1*data.P(i));
        rAMPf(i,j)=filter(b,a,rAMP(i,j));
        rAMP(1+ws:end,j)=rAMPf(1+ws:end,j);
    end
end

data.TOF=TOF;
data.SOS=SOS;
data.rAMP=rAMP;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                          Visualisation                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Visualisation globale des enveloppes des signaux
figure
imagesc(abs(SigHf/max(max(SigHf))))
title('Visualisation compl�te des signaux US')


% % Trac� des zeros crossing et des max sur le premier signal
figure;
plot(S,'b');hold on;plot(tfE,S(tfE),'xr','MarkerSize',10)
plot(tm1(1,:),Am1(1,:),'v')
plot(tm2(1,:),Am2(1,:),'v')
plot(tz1(1,:),Az1(1,:),'s')
plot(tz2(1,:),Az2(1,:),'s')
xlabel('Temps d un tir en nb de point')
ylabel('Amplitude UA')

% %  Trac� des zc et des max sur tous les signaux trait�s
figure
hold on
m=max(max(Sigf));
for i=1:nS
    plot(Sigf(i,:)/m-i)
    plot(SigHf(i,:)/m-i)
    plot(tm1(i,:),Am1(i,:)/m-i,'r+')
    plot(tm2(i,:),Am2(i,:)/m-i,'r+')
    plot(tz1(i,:),Az1(i,:)/m-i,'b+')
    plot(tz2(i,:),Az2(i,:)/m-i,'b+')
end
xlabel('Temps d un tir US en nb de point')
ylabel('Amplitude UA et num�ro du tir')

% % Trac� des spectres
figure;
imagesc(f,tf,abs(spectreS1));
xlim([0 1E7])
title('Spectre echo1')

figure;
imagesc(f,tf,abs(spectreS2));
title('Spectre echo2')
xlim([0 1E7])

figure;
subplot(2,1,1)
plot(tf,fpic1(I,1),'.g',tf,fpic2(I,1),'.r')
xlabel('Temps de compression en sec')
ylabel('fr�quence en Hz')
subplot(2,1,2)
plot(Pf,fpic1(I,1),'.g',Pf,fpic2(I,1),'.r')
xlabel('Hauteur TAXT du mobile en mm')
ylabel('fr�quence en Hz')

figure;
subplot(2,1,1)
plot(tf,apic1(I,1),'.g',tf,apic2(I,1),'.r')
xlabel('Temps de compression en sec')
ylabel({'Amplitude de ','la fr�quence max UA'})
subplot(2,1,2)
plot(Pf,apic1(I,1),'.g',Pf,apic2(I,1),'.r')
xlabel('Hauteur TAXT du mobile en mm')
ylabel({'Amplitude de ','la fr�quence max UA'})



% % Comparaison donn�es TAXT et donn�es US
figure
hold on
plot(tf,data.P./max(data.P))
plot(tf,data.tz2(:,1)./max(data.tz2(:,1)))
plot(tf(tf==0),data.P(tf==0)./max(data.P),'x')
plot(tf(tf==0),data.tz2(tf==0,1)./max(data.tz2(:,1)),'+')
xlabel('Temps de compression en sec')
ylabel('Distance/Tps de vol normalis� UA')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end