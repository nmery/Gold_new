load('ddh_Gold.dat');
load('ddh_sim.dat');
load('Real_Gold.dat');
load('real_sim.dat');
load('ddh_Gold.dat');
load('BM_uc.dat');
load('BM_uc_sim.dat');

x=ddh_Gold;            % ddh dataset
x_dec=ddh_Gold_dec;    % ddh declustered dataset
G=Real_Gold(:,1:3);    % Grid to perform the predictions
R=Real_Gold;           % Real value (interpolated by IDW of RC)
G_UC=BM_uc;            % Panel model for UC

block=[20 10 5];
nd=[3 3 3];
ival=0;
nk_ok=1;
nk_ck=3;
nk_uc=50;
nd_uc=[10 10 5];
rad=99999;
ntok=1;
avg=mean(x(:,4));
smu=[20 10 5];
panneau=[200 100 50];
vc=[0 0.05:0.01:1.5];
dec=0;                 % dec=0 non-declustered histogram for UC / dec=1 declustered histogram for UC

%cas=1 Base case (nk_ok=1 and nk_ck=3)
%cas=2 Same neighborhood for OK and CK (nk=2)
%cas=3 30% increase nug/sill ratio
%cas=4 30% increase range
%cas=5 30% reduce range
%cas=6 Declustered histogram for UC
%cas=7 Reality obtained by turning bands

vcas=[1 7];

for i=1:length(vcas)
    cas=vcas(i)
    
    switch cas
        case 1
            model=[ 1 1; 4 20; 4 80];
            c=[ 0.1; 1.1; 0.1];
            
        case 2
            model=[ 1 1; 4 20; 4 80];
            c=[ 0.1; 1.1; 0.1];
            nk_ok=2;
            nk_ck=2;
            
        case 3
            model=[ 1 1; 4 20; 4 80];
            c=[ 0.3; 0.8 ; 0.1];
            
        case 4
            model=[ 1 1; 4 26; 4 104];
            c=[ 0.1; 1.1 ; 0.1];
            
        case 5
            model=[ 1 1; 4 14; 4 56];
            c=[ 0.1; 1.1; 0.1];
            
        case 6
            model=[ 1 1; 4 20; 4 80];
            c=[ 0.1; 1.1; 0.1];
            dec=1;
            
        case 7
            model=[ 1 1; 4 20; 4 80];
            c=[ 0.1; 1.1; 0.1];
            x=ddh_sim;
            G=real_sim(:,1:3);
            R=real_sim;
            G_UC=BM_uc_sim;
            block=[5 5 5];
            smu=[5 5 5];         
    end
    
    [stat,x0s_ok,x0s_ck,ton_uc,mean_uc]=okckuc(x,x_dec,G,R,model,c,block,nd,nd_uc,ival,nk_ok,nk_ck,nk_uc,rad,ntok,avg,smu,panneau,vc,G_UC,cas,dec);
    
    s_OK=latex(vpa(sym(stat.OK),5))
    s_CK=latex(vpa(sym(stat.CK),5))
    mean_uc(1,1)
    v_sum(i,:)=stat.summary;
       
end