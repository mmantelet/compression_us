function [tm1n, tz1n]=suivi_sig(S, tm1, tz1, fen)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Description                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Programme permettant de d�tecter les points caract�ristiques d'un signal
% (zeros crossing, maximums) � partir des coordonn�es de ces points
% caract�ristiques provenant du tir US pr�c�dent.
% % % % % % % % % % % % % % % Entr�es % % % % % % % % % % % % % % % % % % %
% S : le tir o� l'on doit d�tecter les points caract�ristiques
% (1 X taille Signal)
% tm1 : coordonn�es des maximums
% (1 X n max)
% tz1 : coordonn�es des zeros crossing
% (1 X n zc)
% fen : taille de la fen�tre o� l'on doit chercher le point
% caract�ristique
% (1 X 1)
% % % % % % % % % % % % % % % Sorties % % % % % % % % % % % % % % % % % % %
% tm1n : nouvelles coordonn�es des maximums.
% tz1n : nouvelles coordonn�es des zeros crossing.
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