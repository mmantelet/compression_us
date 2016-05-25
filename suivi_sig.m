function [tm1n, tz1n]=suivi_sig(S, tm1, tz1, fen)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Description                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Programme permettant de détecter les points caractéristiques d'un signal
% (zeros crossing, maximums) à partir des coordonnées de ces points
% caractéristiques provenant du tir US précédent.
% % % % % % % % % % % % % % % Entrées % % % % % % % % % % % % % % % % % % %
% S : le tir où l'on doit détecter les points caractéristiques
% (1 X taille Signal)
% tm1 : coordonnées des maximums
% (1 X n max)
% tz1 : coordonnées des zeros crossing
% (1 X n zc)
% fen : taille de la fenêtre où l'on doit chercher le point
% caractéristique
% (1 X 1)
% % % % % % % % % % % % % % % Sorties % % % % % % % % % % % % % % % % % % %
% tm1n : nouvelles coordonnées des maximums.
% tz1n : nouvelles coordonnées des zeros crossing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tm1n=NaN(length(tm1),2);
tz1n=NaN(length(tz1),2);

for nP=1:length(tm1)
    vec_t=(round(tm1(nP)-fen/2):round(tm1(nP)+fen/2));
    [~, I]=max(abs(S(vec_t)));
    tm1n(nP,1)=vec_t(I);
    tm1n(nP,2)=S(vec_t(I));
    
    
    if nP<=length(tz1)
        vec_z=(round(tz1(nP)-fen/2):round(tz1(nP)+fen/2));
        [~, Imi]=min(S(vec_z(1:(end-1))).*S(vec_z(2:end)));

        p=(S(vec_z(Imi+1))-S(vec_z(Imi)));
        b=S(vec_z(Imi))-p*vec_z(Imi);

        tz1n(nP,1)=-b/p;
        tz1n(nP,2)=0;
    end
    
end



end