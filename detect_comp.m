function [debut, fin]=detect_comp(Vec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Description                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Programme permettant de détecter le début et la fin d'une compression en
% détectant les changement de pente dans le vecteur.
% % % % % % % % % % % % % % % Entrées % % % % % % % % % % % % % % % % % % %
% Vec : Vecteur (Position ou temps )où l'on doit détecter l'inflexion 
% caractéristique du début de la compression.
% % % % % % % % % % % % % % % Sorties % % % % % % % % % % % % % % % % % % %
% debut : début de la compression détecté
%  fin : fin de la compression détecté
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long
fen=1100;
D_vec=abs(diff(Vec));
D_vecf=fft(D_vec);
D_vecf((round(end/2-fen/2)):(round(end/2+fen/2)))=0;    % filtrage du bruit dans le domaine de Fourier
D_Vec=real(ifft(D_vecf));

seuil=max(D_Vec(1:round(end/2)))/2;


debut=find(D_Vec>seuil,1,'first');
f=find(D_Vec(debut:end)<seuil,1,'first');
fin=f+debut;


end