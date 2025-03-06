clear
dd=0;

date=20250208;

filename1=sprintf('%d_cca_mean_rms_all.csv', date);
filename2=sprintf('%d_cca_std_rmspower_all.csv', date);
filename3=sprintf('%d_cta2_mean_rms_all.csv', date);
filename4=sprintf('%d_cta2_std_rmspower_all.csv', date);
filename10=sprintf('%d_y_yD_cca_cta2_mean_rms.csv', date);
sfname=sprintf('%d_200_y_yD_cca_cta2_mean_rms.mat', date);

Uinf=3.99;
sfreq=50000; % Hz
stime=30; % sec
yc=175.22;
yy=[-4 -3 -2 -1 -0.6 -0.4:0.2:1.6 2:0.4:4.0];
yyy=[0];

% P1=[-294.953683 4453.75079 -25161.57506 63053.13974 -59150.11908];
C=1.76226675996221;
D=-10.561264961291;
E=25.5032945300191;
F=-27.5008560512424;
G=11.7014315799771;


P2=[3.643823059	-1.938934606];%較正定数

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for bb=1:22%計測点の数
close all;

yyy_new=yc - 0.3*yy(1,bb);
yyy(bb, 1) = yyy_new;
matobj=matfile(sprintf('%dccacta2_%03d.mat',date, bb));
whos(matobj);
y=matobj.EE(:,1);%ch1
y1=matobj.EE1(:,1);%ch2

sfname5 = sprintf('%d_cca%03d.fig', date, bb);
sfname6 = sprintf('%d_cca%03d.png', date, bb);
sfname7 = sprintf('%d_cta2_%03d.fig', date, bb);
sfname8 = sprintf('%d_cta2_%03d.png', date, bb);

EE_CCA=matobj.EE;%ch1
EE_CTA2=matobj.EE1;%ch2

NN=length(EE_CCA(:,1));
time=[0:NN-1]/sfreq;
NN2=sfreq;
ff=[1:NN2/2]/(NN2)*sfreq;
UU=C*(EE_CCA.^4)+ D*(EE_CCA.^3)+ E*(EE_CCA.^2)+ F*(EE_CCA.^1)+G;%CCA
% U2=UU/Uinf;
U2=UU;

n=0.46;
UU_1=polyval(P2,EE_CTA2.^2).^(1/n); %polyval(P,EE);
% U2_1=UU_1/Uinf;
U2_1=UU_1;



%%%%%%%%%%%%%%%%%%%%%%%流速から実効値の計算%%%%%%%%%%%%%%%%%%%
eff_U1=mean(UU);
eff_U2=std(U2);

eff_U1_1=mean(UU_1);
eff_U2_1=std(U2_1);%これとパワースペクトルと実効値の値がほぼ一緒になる必要がある

cc=0;
for kk=1:1000
    idx=[1:NN2]+(kk-1)*sfreq;
    [idx(end) NN];
    if(idx(end)>NN)
        break;
    else
        cc=cc+1;
        aa=U2(idx);
        aaf=fft(aa)/NN2;
        aafp=aaf.*conj(aaf);
        aafp2=2*aafp(2:NN2/2+1);
        if(kk==1)
            aafpSUM=aafp2;
        else
            aafpSUM=aafpSUM+aafp2;
        end
    end
end
aafpAVE=aafpSUM/cc;

cc=0;
for kk=1:1000
    idx=[1:NN2]+(kk-1)*sfreq;
    [idx(end) NN];
    if(idx(end)>NN)
        break;
    else
        cc=cc+1;
        aa1=U2_1(idx);
        aaf1=fft(aa1)/NN2;
        aafp1=aaf1.*conj(aaf1);
        aafp2_1=2*aafp1(2:NN2/2+1);
        if(kk==1)
            aafpSUM1=aafp2_1;
        else
            aafpSUM1=aafpSUM1+aafp2_1;
        end
    end
end
aafpAVE1=aafpSUM1/cc;

figure(2)
%subplot(1,2,1)
loglog(ff,aafpAVE,'-')
xlabel('{\it f} (Hz)')
ylabel('{\it Power}')
grid on
title(sprintf('CCA y=%.2f', yyy))
set(gca,'Fontname','Times','Fontsize',14)
% saveas(gcf, sfname5);
% saveas(gcf, sfname6);

figure(3)
%subplot(1,2,1)
loglog(ff,aafpAVE1,'-')
xlabel('{\it f} (Hz)')
ylabel('{\it Power}')
grid on
title(sprintf('CTA2 y=%.2f', yyy))
set(gca,'Fontname','Times','Fontsize',14)
% saveas(gcf, sfname7);
% saveas(gcf, sfname8);

%%%%%%%%%%%%パワースペクトルから実効値の計算%%%%%%%%%%%%%%

aafpAVESUM=sum(aafpAVE);
eff_spc=aafpAVESUM^(0.5)

aafpAVESUM1=sum(aafpAVE1);
eff_spc_1=aafpAVESUM1^(0.5)

Urms_all(bb, 1)=eff_U2;
Umean_all(bb, 1)=eff_U1;
Urms_power_cca(bb, 1)=eff_spc;

Urms_all1(bb, 1)=eff_U2_1;
Umean_all1(bb, 1)=eff_U1_1;
Urms_power_cta(bb, 1)=eff_spc_1;

end

combined_data=[Umean_all Urms_all];
combined_data1=[Urms_all Urms_power_cca];
writematrix(combined_data,filename1);
writematrix(combined_data1,filename2);

combined_data2=[Umean_all1 Urms_all1];
combined_data3=[Urms_all1 Urms_power_cta];
writematrix(combined_data2,filename3);
writematrix(combined_data3,filename4);

combined_data_all=[yyy yy' Umean_all Urms_all Umean_all1 Urms_all1];
writematrix(combined_data_all,filename10);
yD=yy';
save(sfname, 'yyy', 'yD', 'Umean_all', 'Urms_all', 'Umean_all1', 'Urms_all1');
disp('END')
