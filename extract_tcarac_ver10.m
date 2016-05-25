function [echo_win, tm, tz, t1z, t1m, peaks,locs,widths]=extract_tcarac_ver10(data,seuil1,seuil2,T_fen, nP, coeff_nbp,fig)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%                            Description                                  %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Fonction permettant de récupérer les premiers passages à zeros, 
% les amplitudes et les fenêtres de chaque écho de chaque échos du tir.
% % % % % % % % % % % % % % % % Entrées % % % % % % % % % % % % % % % % % %
% data = 1 tir ultrasonore
% fig=1 : affichage des figures, fig=0 : pas d'affichage de figure.
% seuil1= seuil de détection des échos selon la hauteur
% seuil2= seuil de détection suivant la largeur
% T_fen = taille de la fenêtre/2 pour la detection des pics de Hilbert
% coeff_nbp = coeff multiplicateur du nombre de points temporels
% nP : nombre de max dans un écho à récupérer
% % % % % % % % % % % % % % % Sorties % % % % % % % % % % % % % % % % % % %
% echo_win = coordonnées des fenêtres temproelles de chaque échos (
% echo_win(:,1)=début des échos et echo_win(:,2)=fin des échos )
% tm = structure comportant les positions et les amplitudes des maximums
% ex :      tm.echo(n echo).max(nPic,1(position) 2(amplitude))
% tz = structure comportant les positions et les amplitudes des zeros
% crossing
% ex :      tz.echo(n echo).zc(nPic,1(position) 2(amplitude))
% t1m = valeur du premier maximum de chaque échos (dim 1 X numéro de
% l'écho)
% t1z = valeur du premier passage à zeros de chaque échos
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%%
% % % % % % % % % % % Prétraitement du tir % % % % % % % % % % % % % % % % 
% nP=5;%nombre de max à recueillir
data=data-mean(data);

%%
% % % % % % % % % Détection de contour Hilbert % % % % % % % % % % % % % % % 

Hil=abs(hilbert(data));
Adata=abs(data);



[peaks,locs,widths,~]=findpeaks(Hil,1:length(Hil),'MinPeakHeight',seuil1*max(data),'MinPeakProminence',0.01,'MinPeakDistance',seuil2,'Annotate','extents');


start_pulse=locs;%-round(seuil2*coeff_nbp);

end_pulse=locs;%+round((seuil2+70)*coeff_nbp);

% T_fen=50;%taille de la fenêtre d'echo


for ne=1:length(start_pulse)
    if (start_pulse(ne)-round(T_fen*coeff_nbp))>=1
        d_hil_s=diff(Hil((start_pulse(ne)-round(T_fen*coeff_nbp)):start_pulse(ne)));
        start_pulse(ne)=start_pulse(ne)-round(T_fen*coeff_nbp);
    else
        start_pulse(ne)=1;
        d_hil_s=diff(Hil(1:start_pulse(ne)));
    end
    if (round((T_fen+80)*coeff_nbp+end_pulse(ne)))<length(Hil)
       d_hil_e=diff(Hil((end_pulse(ne)):(round(T_fen*coeff_nbp)+end_pulse(ne))));
       end_pulse(ne)=round((T_fen+80)*coeff_nbp+end_pulse(ne));
    else
        end_pulse(ne)=length(data);
       d_hil_e=diff(Hil((end_pulse(ne)):end));
    end
end



echo_win=[start_pulse', end_pulse'];

%%
% % % % Détection des max et des passages à zeros de chaque fenêtre % % % % 

for n_start=1:length(start_pulse)

    vec_t=(start_pulse(n_start)):(end_pulse(n_start));

    max_pulse=max(Adata(vec_t));
    [~,locsm]=findpeaks(Adata(vec_t),'MinPeakHeight',0.025*max_pulse,'MinPeakDistance',0.1*T_fen,'annotate','peaks');
    tm.echo(n_start).max(:,1)=vec_t(locsm(1:nP));
    tm.echo(n_start).max(:,2)=data(vec_t(locsm(1:nP)));


    if ~isempty(locsm)
            vec_z=(tm.echo(n_start).max(1,1)):(tm.echo(n_start).max(end,1));%(end_pulse(n_start));
        cz=0;
        for z=1:length(vec_z)
            if data(vec_z(z)-1)*data(vec_z(z))<0
                cz=cz+1;
                p=(data(vec_z(z))-data(vec_z(z)-1))/(vec_z(z)-(vec_z(z)-1));
                b=data(vec_z(z))-p*vec_z(z);
                
                tz.echo(n_start).zc(cz,1)=(-b/p);
                tz.echo(n_start).zc(cz,2)=0;
            end
        end
    end
end

t1m(:,1)=tm.echo(1).max(1,1);
t1m(:,2)=tm.echo(1).max(1,2);

t1z(:,1)=tz.echo(1).zc(1,1);
t1z(:,2)=tz.echo(1).zc(1,2);


%%
% % % % % % % % % % % % % % Visualisation % % % % % % % % % % % % 

if fig==1
    figure;
    plot((data),'g');
    hold on
    plot(abs(data),'--g')
    hold on
    plot(abs(hilbert(data)),'-xr')
    hold off
    title('Visualisation signal')
    legend('data','abs(data)','abs(hilbert(data))')


    figure;
    plot(data,'-b')
    hold on
    plot(start_pulse,data(start_pulse),'dg','MarkerSize',10)
    plot(end_pulse,data(end_pulse),'dr','MarkerSize',10)
    for n_start=1:length(tm.echo)
        plot(tm.echo(n_start).max(:,1),tm.echo(n_start).max(:,2),'xr','MarkerSize',10)
        plot(tz.echo(n_start).zc(:,1),tz.echo(n_start).zc(:,2),'xg','MarkerSize',10)

    end
    plot(1:length(data),ones(length(data),1)*seuil1,'--r')
    plot(start_pulse'*ones(1,101),(-1:0.02:1),'--c')
    plot(end_pulse'*ones(1,101),(-1:0.02:1),'--m')
    hold off
end


end