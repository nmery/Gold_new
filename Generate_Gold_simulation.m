load('mayo_data_for_denis.mat');

cap=30;  %capping value

% Capping

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

% Defining domain of RC

xmin=floor(min(rc(:,1)));
xmax=round(max(rc(:,1)));
dx=5;
ymin=floor(min(rc(:,2)));
ymax=round(max(rc(:,2)));
dy=5;
zmin=floor(min(rc(:,3)));
zmax=round(max(rc(:,3)));
dz=5;

BM_rc=grille3(xmin,xmax,dx,ymin,ymax,dy,zmin,zmax,dz);

x0=BM_rc;
x=rc;
a=0.000000000000000000000000001;
dmax=6.3;
nmax=240;
[x0s_rc,somme,~,~]=invd_new_3d(x,x0,a,dmax,nmax);

x0s_rc(any(isnan(x0s_rc),2),:)=[];


% Obtain anamorphosis of rc

y=anamor(rc(:,4));

% Apply TB

x0=x0s_rc(:,1:3);
x=[rc(:,1:3),y(:,2)];
model=[1 1 ; 2 15];
c=[0.57;0.36];
seed=12345678;
nbsimul=1;

datasim=tourband_cond_large(x,x0,model,c,seed,nbsimul);

% Post processing TB

yz=[rc(:,1:3),y];
[z]=anamorinv(yz,datasim(:,4));
real=[x0,z(:,2)];

save('real_sim.dat','real','-ascii')
stats_Real=[length(real(:,4)),mean(real(:,4)),var(real(:,4)),(std(real(:,4))/mean(real(:,4)))];
s_Real=latex(vpa(sym(stats_Real),5))

% Define ddh dataset

x1=real;
x01=ddh(:,1:3);
a=1;
dmax=6.5;
nmax=1;

[x0s_d,~,~]=invd_new_3d(x1,x01,a,dmax,nmax);

x0s_d(any(isnan(x0s_d),2),:)=[];
idx=ismember(ddh(:,1:3),x0s_d(:,1:3),'rows');
ddh_dat=idx.*ddh;
ddh_dataset= ddh_dat(any(ddh_dat,2),:);

stats_ddh_dataset=[length(ddh_dataset(:,4)),mean(ddh_dataset(:,4)),var(ddh_dataset(:,4)),(std(ddh_dataset(:,4))/mean(ddh_dataset(:,4)))];
s_ddh_dataset=latex(vpa(sym(stats_ddh_dataset),5))
save('ddh_sim.dat','ddh_dataset','-ascii')

% Panel model for UC

xmin=floor(min(real(:,1)));
xmax=round(max(real(:,1)));
dx=200;
ymin=floor(min(real(:,2)));
ymax=round(max(real(:,2)));
dy=100;
zmin=floor(min(real(:,3)));
zmax=round(max(real(:,3)));
dz=50;

xmin=xmin-2.5+100;
xmax=xmax+2.5-100;
ymin=ymin-2.5+50;
ymax=ymax+2.5-50;
zmin=zmin-2.5+25;
zmax=zmax+2.5-25;

BM_uc=grille3(xmin,xmax,dx,ymin,ymax,dy,zmin,zmax,dz);

save('BM_uc_sim.dat','BM_uc','-ascii')