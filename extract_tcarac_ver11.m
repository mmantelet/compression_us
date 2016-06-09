function [echo_win,Techo,tm, tz]=extract_tcarac_ver11(S,Sref,T_fen,nP,fig)
%%
warning('off')
% Sinit=A(750:end,1);
F=1E6;                  %Fréquence du pulse
% T=1/F*1*(100E6);        %Période temporelle du signal
% nP=5;                   %Nombre de pics à détecter
% nE=2;                   %Nombre d'échos analysés
% seuil1=0.4;             %hauteur minimal de détection des pics sur la transformée de Hilbert. Paramètre compris entre 0 et 1
% seuil2=200;             %écart minimal entre 2 pics sur la transformée de Hilbert
% coeff_nbp=1;            %coefficient multiplicateur du nombre de point du signal (en cas d'interpolation coeff_nbp >1)
echo_fen=[100, 200];    %coordonnée de la fenêtre d'un echo par rapport à zc1 et zc4. ex de coord [zc1-echo_fen(1), zc4+echo_fen(2)]
n_inter=10;
% T_fen=3*T;
% [echo_win, tm, tz, t1z, t1m, peaks,locs,widths]=extract_tcarac_ver10(Sinit,seuil1,seuil2,3*T, nP, coeff_nbp, 1);
% % Sref=Sinit(tz.echo(2).zc(2,1):tz.echo(2).zc(end,1),1);
% load('Srefeau.mat')
% Sref=Srefeau;


    [fctcorr, retard]=xcorr((S),(Sref));
    Scorr=(fctcorr)./max(fctcorr);

    [~,l,widths,~] = findpeaks(abs(hilbert(Scorr)),'MinPeakProminence',0.05);
    
    
    f_inf1=round(l(1)-widths(1)/2);
    f_sup1=round(l(1)+widths(1)/2);
    [~, I]=max(Scorr(f_inf1:f_sup1));
    milS1=retard(I(1)+l(1));
    
    xinterp1(1,:)=retard((I(1)+f_inf1-n_inter):(I(1)+f_inf1+n_inter));
    yinterp1(1,:)=Scorr((I(1)+f_inf1-n_inter):(I(1)+f_inf1+n_inter));
    pp1=polyfit(xinterp1,yinterp1,2);
    t1_intercorr=-pp1(2)./(2*pp1(1));


    f_inf2=round(l(2)-widths(2)/2);
    f_sup2=round(l(2)+widths(2)/2);
    [~, I2]=max(Scorr(f_inf2:f_sup2));
    milS2=retard(I2(1)+l(2));
    
    xinterp2(1,:)=retard((I2(1)+f_inf2-n_inter):(I2(1)+f_inf2+n_inter));
    yinterp2(1,:)=Scorr((I2(1)+f_inf2-n_inter):(I2(1)+f_inf2+n_inter));
    pp2=polyfit(xinterp2,yinterp2,2);
    t2_intercorr=-pp2(2)./(2*pp2(1));
    
%     TOF_intercorr(i)=t2_intercorr-t1_intercorr;
    
    
    Techo(1,1)=t1_intercorr;
    Techo(1,2)=t2_intercorr;

    
    f_S_inf1=round(milS1(1)-length(Sref)/2);
    f_S_sup1=round(milS1(1)+length(Sref)/2);
    f_S_inf2=round(milS2(1)-length(Sref)/2);
    f_S_sup2=round(milS2(1)+length(Sref)/2);
   
    
    
    
%     [M3, I3]=max(Scorr(round(l(3)-widths(2)/2):round(l(3)+widths(2)/2)));
%     milS3=retard(I3(1)+l(3));
%     plot(round(milS3(1)-length(Sref/2)):round(milS3(1)+length(Sref/2)),S(round(milS3(1)-length(Sref/2)):round(milS3(1)+length(Sref/2))),'m')

%% Détection max et passage à zeros

start_pulse=[f_S_inf1, f_S_inf2]-echo_fen(1);
end_pulse=[f_S_sup1, f_S_sup2]+echo_fen(2);
echo_win(:,1)=start_pulse';
echo_win(:,2)=end_pulse';
Adata=abs(S);

for n_start=1:length(start_pulse)

    vec_t=(start_pulse(n_start)):(end_pulse(n_start));

    max_pulse=max(Adata(vec_t));
    [~,locsm]=findpeaks(Adata(vec_t),'MinPeakHeight',0.05*max_pulse,'MinPeakDistance',0.1*T_fen,'annotate','peaks');
    tm.echo(n_start).max(:,1)=vec_t(locsm(1:nP));
    tm.echo(n_start).max(:,2)=S(vec_t(locsm(1:nP)));


    if ~isempty(locsm)
        cz=1;cano=0;
        for nmax=1:(length(tm.echo(n_start).max(:,1))-1);
                
                vec_z=(tm.echo(n_start).max(nmax,1)):(tm.echo(n_start).max(nmax+1,1));%(end_pulse(n_start));
                
                czi=cz;
            
                for z=1:length(vec_z)
                    if S(vec_z(z)-1)*S(vec_z(z))<0 || S(vec_z(z))==0

                        p=(S(vec_z(z))-S(vec_z(z)-1))/(vec_z(z)-(vec_z(z)-1));
                        b=S(vec_z(z))-p*vec_z(z);

                        tz.echo(n_start).zc(cz,1)=(-b/p);
                        tz.echo(n_start).zc(cz,2)=0;
                        cz=cz+1;
                        moni(1,1)=czi;moni(1,2)=cz;
                    end
                end
            
                if abs(czi-cz)>1;
                        cano=cano+1;
                        pano(cano,1)=nmax;
                        pano(cano,2)=cz;
                        pano(cano,3)=czi;
                        pano(cano,4)=length(tz.echo(n_start).zc(:,1));
                        tz.echo(n_start).zc(czi:end,:)=[];
                        tz.echo(n_start).zc(czi,1)=NaN(1);
                        tz.echo(n_start).zc(czi,2)=NaN(1);
                        cz=czi+1;
                elseif abs(czi-cz)==0
                        tz.echo(n_start).zc(cz,1)=NaN(1);
                        tz.echo(n_start).zc(cz,2)=NaN(1);
                        cz=cz+1;
                end

        end
        
        if (length(tz.echo(n_start).zc(:,1))-cz)~=0
            tz.echo(n_start).zc(cz+1:length(tz.echo(n_start).zc(:,1)),:)=[];
        end
    end
end

t1m(:,1)=tm.echo(1).max(1,1);
t1m(:,2)=tm.echo(2).max(1,2);

t1z(:,1)=tz.echo(1).zc(1,1);
t1z(:,2)=tz.echo(2).zc(1,2);
    
%% Visualisation

if fig==1
   figure
   subplot(2,1,1)
   plot(S);
   hold on
   plot(start_pulse(1):end_pulse(1),S(start_pulse(1):end_pulse(1)),'.r')
   plot(start_pulse(2):end_pulse(2),S(start_pulse(2):end_pulse(2)),'.c')
   plot(tm.echo(1).max(:,1),tm.echo(1).max(:,2),'or','MarkerFaceColor','k','MarkerSize',5)
   plot(tm.echo(2).max(:,1),tm.echo(2).max(:,2),'or','MarkerFaceColor','k','MarkerSize',5)
   plot(tz.echo(1).zc(:,1),tz.echo(1).zc(:,2),'square b','MarkerFaceColor','k','MarkerSize',5)
   plot(tz.echo(2).zc(:,1),tz.echo(2).zc(:,2),'square b','MarkerFaceColor','k','MarkerSize',5)

   hold off
    
   subplot(2,1,2)
   findpeaks(abs(hilbert(Scorr)),retard,'MinPeakProminence',0.05)
    xlim([0 4500])
end

end