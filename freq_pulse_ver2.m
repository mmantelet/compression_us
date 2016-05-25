function [freq, A6dB, apic, Dfreq]=freq_pulse_ver2(data, echo_win, coeff_nbp, fig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Description                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Programme permettant de r�cup�rer toutes les donn�es fr�quentielles
% (spectre, bande passante,...)
% % % % % % % % % % % % % % % Entr�es % % % % % % % % % % % % % % % % % % %
% data : tir � analyser
% (1 X taille Signal)
% echo_win : coordonn�es des fen�tres des �chos � analyser 
% [debut echo, fin echo]
% coeff_nbp : coefficient multiplicateur du nombre de points dans un tir
% (1 X 1)
% fig : 1 affichage des figures, 0 pas d'affichage.
% % % % % % % % % % % % % % % Sorties % % % % % % % % % % % % % % % % % % %
% freq : fr�quence du de la fr�quence maximale
% A6dB : largeur � 6dB (mi-hauteur) de la distribution en fr�quence
% apic : l'amplitude de la fr�quence maximale.
% Dfreq : structure comportant toutes les donn�es fr�quentielles (spectres,
% fr�quence, phase, ...)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning('off')

[li, ~]=size(echo_win);
n_f=1024;%length(data);
for pul=1:li;
    masque=ones(size(data(echo_win(pul,1):echo_win(pul,2))));%tukeywin(length(data(echo_win(pul,1):echo_win(pul,2))),0.2);
    data_f=data(echo_win(pul,1):echo_win(pul,2)).*masque;
    fftdataf=(fft(data_f,2028));

    fe=(10^8)*coeff_nbp;
    f=(1:length(fftdataf))*fe./length(fftdataf);
    t=(1:length(data_f))./fe;
    
    [apic, fpic, A6dB, ~]=findpeaks(abs(fftdataf),1:length(f),'MinPeakHeight',0.5*max(abs(fftdataf)),'MinPeakProminence',1,'Annotate','extents');%'Threshold',0.5);%max(fftdataf)/2);

    
    if ~isempty(fpic)
        freq(1,pul)=f(round(fpic(1)));
    else
        freq(1,pul)=NaN(1,1);
    end
    
    ang=unwrap(angle(fftdataf));
    
    
    Dfreq.echo(pul).echowin=echo_win(pul,1):echo_win(pul,2);
    Dfreq.echo(pul).S=data_f;
    Dfreq.echo(pul).t=t;
    Dfreq.echo(pul).fftS=fftdataf;
    Dfreq.echo(pul).phase=ang;
    Dfreq.echo(pul).f=f;
    Dfreq.echo(pul).fpic=freq;
    Dfreq.echo(pul).apic=apic;
    Dfreq.echo(pul).A6dB=A6dB;
    
end


if fig==1
    figure;
    plot(t,data_f)

    figure;
    plot(f,abs(fftdataf))

    figure;
    plot(f',unwrap(angle(fftdataf)))
end

end