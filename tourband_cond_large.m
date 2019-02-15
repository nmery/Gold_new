function datasim=tourband_cond_large(x,x0,model,c,seed,nbsimul)
% fonction pour r�aliser simulation conditionnelle par bandes tournantes de
% moyenne 0
% syntaxe : datasim=tourband_cond(x,x0,model,c,seed,nbsimul)
%
% x: matrice n x 4 des donn�es conditionnantes
% x0: matrice nx0 x 3 des points ou r�aliser les simulations
% model et c comme cokri mais seulement des modeles isotropes
% seed : germe al�atoire
% nbsimul: nombre de r�alisations
%
% ATTENTION: le programme suppose que les donn�es en x suivent une loi
% normale de moyenne 0 (ne fait pas la transformation) et que celle-ci est compatible
% avec le modele de variogramme fourni.
% Chaque composante de model est isotrope.
%
% Auteur: D. Marcotte novembre 2004



% g�n�rer les r�alisations n.c. aux points �chantillons et aux points de grille
nx0=size(x0,1);
ns=size(x,1);
nlim=500;     % nombre maximum de points ou blocs � consid�rer simultan�ment

x0=tourband([x0;x(:,1:3)],model,c,seed,nbsimul); 

xsim=x0(nx0+1:end,:);
x0=x0(1:nx0,:);

% effectuer le post-conditionnment par krigeage
if ns>nlim | nx0>nlim
      
    nsteps=ceil(nx0/nlim);       % on d�coupe en nlim blocs
    for i=1:nsteps
        ideb=(i-1)*nlim+1;
        ifin=min(i*nlim,nx0);
        xx0=x0(ideb:ifin,:); 
        idx=knnsearch(x(:,1:3),xx0(:,1:3),'K',10);
        idx=unique(idx(:));          % idx continet les points utilis�s pour estimer la s�rie de blocs actuels
        [i,length(idx)];
        datasim(ideb:ifin,:)=postcond(x(idx,:),xsim(idx,:),xx0,model,c,nbsimul);
    end
else
   datasim=postcond(x,xsim,x0,model,c,nbsimul);
end



function datasim=postcond(x,xsim,x0,model,c,nbsimul)
% fonction pour r�aliser le post-conditionnement
% suppose x, xsim et x0 de moyenne 0
nx0=size(x0,1);
np=500; % on conditionne les points par groupe de 500
% calculer matrice de covariance
k=covardm(x(:,1:3),x(:,1:3),model,c);
ki=inv(k); clear k;

b0=ki*x(:,end);
datasim=x0;
bi=ki*xsim(:,4:end);
for j=1:np:nx0;
   ifin=min(j-1+np,nx0);
    k0=covardm(x(:,1:3),x0(j:ifin,1:3),model,c);
    datasim(j:ifin,4:end)=((b0*ones(1,nbsimul)-bi)'*k0)'+x0(j:ifin,4:end);
end





