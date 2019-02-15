function datasim=tourband(x0,model,c,seed,nbsimul)
% fonction pour réaliser simulation non-conditionnelle par bandes
% tournantes (de moyenne 0)
% syntaxe : datasim=tourband(x0,model,c,seed,nbsimul)
%
% x0: matrice nx0 x 3 des points ou réaliser les simulations
% model et c comme cokri mais seulement des modeles isotropes
% seed : germe aléatoire
% nbsimul: nombre de réalisations
%
% ATTENTION: Chaque composante de model est isotrope.
%
% Auteur: D. Marcotte novembre 2004
%

nl=500; % nombre de lignes
nv=500; % nombre de valeurs sur chaque ligne
randn('state',seed);
rand('state',seed);
nx0=size(x0,1);

x0c=(max(x0)+min(x0))/2;
dmax=(max(x0)-min(x0));dmax=sqrt(sum(dmax.^2)); % distance max rencontrée
dx=dmax/nv;

x0old=x0;
for i=1:3
    x0(:,i)=(x0(:,i)-x0c(i))/dmax*nv; % la grille est centrée réduite;
end

% générer les directions sur la sphere
ligne=dirsphere_van(nl);

model2=model;model2(:,2)=model2(:,2)/dmax*nv; % ajuster la portée a l'échelle 
id=model2(:,1)==1; % identifier l'effet de pépite
if sum(id)==1;
    zst=randn(nx0,nbsimul)*sqrt(c(id));
else
    zst=zeros(nx0,nbsimul);
end
model2=model2(~id,:);c=c(~id);

if size(model2,1)>0; % si ce n'est pas un effet de pépite pur
    
    k=covarligne(model2,c,nv);
    
    % faire LU sur les lignes
    u=chol(k);
    
    % faire une réalisation a la fois
    
    for ir=1:nbsimul;
        % traiter ligne par ligne
        zs=zeros(nx0,1);
        
        for il=1:nl
            
            %projeter les coordonnées sur les lignes
            ipos=ceil(x0*ligne(il,:)')+nv/2;
            
            % imposer la covariance sur une ligne
            ys=u'*randn(nv,1);
            
            % cumuler les valeurs trouvées sur chaque ligne
            zs=zs+ys(ipos);
            
        end
        
        zs=1/sqrt(nl)*zs;
        zst(:,ir)=zst(:,ir)+zs;
    end
end

datasim=[x0old zst];

function ligne=dirsphere(nl);
nl2=round(2*nl+5*sqrt(nl/2));
x=rand(nl2,3)*2-1;
lx=sqrt(sum(x.^2,2));
id=lx<1;
x=x(id,:)./(lx(id)*ones(1,3));
ligne=x(1:nl,:);

function k=covarligne(model,c,nv);

nm=size(model,1);
x=[1:nv]';
h=abs(x*ones(1,nv)-ones(nv,1)*x');
k=zeros(nv,nv);
for i=1:nm
    a=model(i,2);
    switch model(i,1)
        case 1 % pépite 
            % ne rien faire
        case 2 % exponentiel
            k=k+c(i)*exp(-h/a).*(1-h/a);
        case 3 % gaussien
            k=k+c(i)*exp(-(h/a).^2).*(a^2-2*h.^2)/a^2;
            k=rendreposdef(k);
        case 4 %spherique
            h2=min(h,a)/a;
            k=k+c(i)*(1-3*h2+2*h2.^3);
        case 6 % cubique
            h2=min(h,a)/a;
            k=k+c(i)*(1-21*h2.^2+35*h2.^3-21*h2.^5+6*h2.^7);

    end
end

function ligne=dirsphere_van(nl);
for i=1:nl;
    a2=base(i,2);
    a3=base(i,3);
    u(i,1)=sum(a2./(2.^[1:size(a2,2)]));
    u(i,2)=sum(a3./(3.^[1:size(a3,2)]));
end
ligne=[cos(2*pi*u(:,1)).*sqrt(1-u(:,2).^2),sin(2*pi*u(:,1)).*sqrt(1-u(:,2).^2),u(:,2)];


function y=base(n,b);
% trouver k
kmax=floor(log(n)/log(b));
n2=n;
for k=kmax:-1:0
   y(k+1)=floor(n2/b^k);
   n2=mod(n2,b^k);
end




