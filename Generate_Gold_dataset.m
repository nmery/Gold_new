% Input

load('mayo_data_for_denis.mat');
cap=30;  %capping value

%Capping

for i=1:length(rc)
    if rc(i,4)>=cap
        rc(i,4)=cap;
    else
    end
end

for i=1:length(ddh)
    if ddh(i,4)>=cap
        ddh(i,4)=cap;
    else
    end
end

% 1. Define grid size

xmin=floor(min(rc(:,1)));
xmax=round(max(rc(:,1)));
dx=20;
ymin=floor(min(rc(:,2)));
ymax=round(max(rc(:,2)));
dy=10;
zmin=floor(min(rc(:,3)));
zmax=round(max(rc(:,3)));
dz=5;

BM_rc=grille3(xmin,xmax,dx,ymin,ymax,dy,zmin,zmax,dz);

% 2. Inverse distance RC

x0=BM_rc;
x=rc;
a=0.000000000000000000000000001;
dmax=6.3;
nmax=240;

[x0s_rc,somme,~,~]=invd_new_3d(x,x0,a,dmax,nmax);

x0s_rc(any(isnan(x0s_rc),2),:)=[];
save('Real_gold.dat','x0s_rc','-ascii')

stats_Real=[length(x0s_rc(:,4)),mean(x0s_rc(:,4)),var(x0s_rc(:,4)),(std(x0s_rc(:,4))/mean(x0s_rc(:,4)))];
s_Real=latex(vpa(sym(stats_Real),5))

% 3. Ddh dataset

x1=x0s_rc;
x01=ddh(:,1:3);
a=3;
dmax=10.3;
nmax=1;

[x0s_d,~,~]=invd_new_3d(x1,x01,a,dmax,nmax);

x0s_d(any(isnan(x0s_d),2),:)=[];
idx=ismember(ddh(:,1:3),x0s_d(:,1:3),'rows');
ddh_dat=idx.*ddh;
ddh_dataset= ddh_dat(any(ddh_dat,2),:);
save('ddh_gold.dat','ddh_dataset','-ascii')

stats_ddh_dataset=[length(ddh_dataset(:,4)),mean(ddh_dataset(:,4)),var(ddh_dataset(:,4)),(std(ddh_dataset(:,4))/mean(ddh_dataset(:,4)))];
s_ddh_dataset=latex(vpa(sym(stats_ddh_dataset),5))

% 4. Declustered histogram ddh

xmin=floor(min(x0s_rc(:,1)));
xmax=round(max(x0s_rc(:,1)));
dx=3;
ymin=floor(min(x0s_rc(:,2)));
ymax=round(max(x0s_rc(:,2)));
dy=3;
zmin=floor(min(x0s_rc(:,3)));
zmax=round(max(x0s_rc(:,3)));
dz=3;

BM_rc_fine=grille3(xmin,xmax,dx,ymin,ymax,dy,zmin,zmax,dz);

x1=ddh_dataset;
x01=BM_rc_fine;
a=1;
dmax=300;
nmax=1;

[x0s_dec,~,~]=invd_new_3d(x1,x01,a,dmax,nmax);

save('ddh_gold_dec.dat','x0s_dec','-ascii')

% 5. Panel model

dx=200;
dy=100;
dz=50;

xmin=xmin-10+100;
xmax=xmax+10-100;
ymin=ymin-5+50;
ymax=ymax+5-50;
zmin=zmin-2.5+25;
zmax=zmax+2.5-25;

BM_uc=grille3(xmin,xmax,dx,ymin,ymax,dy,zmin,zmax,dz);
save('BM_uc.dat','BM_uc','-ascii')
