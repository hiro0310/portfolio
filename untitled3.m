clear
close all

date=20250205;
t = 60;

filename1=sprintf('%d_cca_mean_rms_max%ds.csv', date, t);
filename2=sprintf('%d_cca_std_rmspower_max%ds.csv', date, t);
filename3=sprintf('%d_cta2_mean_rms_max%ds.csv', date, t);
filename4=sprintf('%d_cta2_std_rmspower_max%ds.csv', date, t);
filenamemax=sprintf('%d_yrmsmax_cca_cta2_mean_rms%ds.csv', date, t);

Uinf=4.005769323;
sfreq=50000; % Hz
stime=t; % sec

yrmsmax=175.52;
yyy=yrmsmax;
%P1=[28.13505375,-276.1672258,1026.990679,-1705.958579,1064.857895];
C=1.45880722333814;
D=-8.2080264406535;
E=18.8716454156307;
F=-18.6449417719788;
G=7.41068069106366;

P2=[3.568998225	-1.76996656];%較正定数
% close all;

matobj=matfile(sprintf('%dccacta2rmsmax%ds_001.mat', date, t));
whos(matobj);
y=matobj.EE(:,1);%ch1
y1=matobj.EE1(:,1);%ch2

sfname5 = sprintf('%d_ccamax%ds.fig', date, t);
sfname6 = sprintf('%d_ccamax%ds.png', date, t);
sfname7 = sprintf('%d_cta2max%%ds.fig', date, t);
sfname8 = sprintf('%d_cta2max%ds.png', date, t);
sfname9 = sprintf('%d_ccacta2max%ds.fig', date, t);
sfname10 = sprintf('%d_ccacta2max%ds.png', date, t);

EE_CCA=matobj.EE;%ch1
EE_CTA=matobj.EE1;%ch2

NN=length(EE_CCA(:,1));
time=[0:NN-1]/sfreq;
NN2=sfreq;
ff=[1:NN2/2]/(NN2)*sfreq;

UU_CCA=C*(EE_CCA.^4)+ D*(EE_CCA.^3)+ E*(EE_CCA.^2)+ F*(EE_CCA.^1)+G;%CCA
U2=UU_CCA/Uinf;

n=0.46;
UU_1=polyval(P2,EE_CTA.^2).^(1/n); %polyval(P,EE);
U2_1=UU_1/Uinf;


%%%%%%%%%%%%%%%%%%%%%%%流速から実効値の計算%%%%%%%%%%%%%%%%%%%
eff_U1=mean(UU_CCA);
eff_U2=std(U2);%不偏標準偏差してる調査実験で得たサンプルデータから母集団のばらつきを推定するときに使用します。

eff_EE=mean(EE_CCA);

eff_U1_1=mean(UU_1);
eff_U2_1=std(U2_1);%これとパワースペクトルと実効値の値がほぼ一緒になる必要がある

%UU2_1のパワースペクトルの平均を計算
% U2_1 を NN2 点ごとに分割（のループ）  
% 各区間ごとに FFT（高速フーリエ変換）を実行
% フーリエ変換後の信号の大きさを計算（パワースペクトル）
% すべての区間のパワースペクトルを合計
% 最後に平均を取ることで、全体の平均パワースペクトル密度 aafpAVE1 を求める

cc=0;
for kk=1:1000
    idx=[1:NN2]+(kk-1)*(sfreq/100);%データ区間の取得sfreqごと
    [idx(end) NN];
    if(idx(end)>NN)%データ範囲が NN を超えたら終了
        break;
    else
        cc=cc+1;
        %        aa=EE(idx,1);
        aa=U2(idx);
        aaf=fft(aa)/NN2;%FFTして正規化
        aafp=aaf.*conj(aaf);
        aafp2=2*aafp(2:NN2/2+1);%パワースペクトルの対称性を考慮　二倍
        if(kk==1)
            aafpSUM=aafp2;
        else
            aafpSUM=aafpSUM+aafp2;%すべてのブロックのパワースペクトルを足し合わせていく
        end
    end
end
aafpAVE=aafpSUM/cc;%すべてのブロックのパワースペクトルを平均ccは30s

cc=0;
for kk=1:1000
    idx=[1:NN2]+(kk-1)*(sfreq/100);
    [idx(end) NN];
    if(idx(end)>NN)
        break;
    else
        cc=cc+1;
        %        aa=EE(idx,1);
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
xlabel('{\it f}(Hz)')
ylabel('{\it Power}')
grid on
title(sprintf('CCA y=%.2f', yrmsmax));
set(gca,'Fontname','Times','Fontsize',14)
saveas(gcf, sfname5);
saveas(gcf, sfname6);

figure(3)
%subplot(1,2,1)
loglog(ff,aafpAVE1,'-')
xlabel('{\it f}(Hz)')
ylabel('{\it Power}')
grid on
title(sprintf('CTA y=%.2f', yrmsmax));
set(gca,'Fontname','Times','Fontsize',14)
saveas(gcf, sfname7);
saveas(gcf, sfname8);


%%%%%%%%%%%%パワースペクトルから実効値の計算%%%%%%%%%%%%%%

aafpAVESUM=sum(aafpAVE);
eff_spc=aafpAVESUM^(0.5)

aafpAVESUM1=sum(aafpAVE1);
eff_spc_1=aafpAVESUM1^(0.5)

bb=0;
Urms_all(bb+1, 1)=eff_U2;
Umean_all(bb+1, 1)=eff_U1;
EEmean_all(bb+1, 1)=eff_EE;
Urms_power_cca(bb+1, 1)=eff_spc;

Urms_all1(bb+1, 1)=eff_U2_1;
Umean_all1(bb+1, 1)=eff_U1_1;
Urms_power_cta(bb+1, 1)=eff_spc_1;


combined_data=[Umean_all Urms_all];
combined_data1=[Urms_all Urms_power_cca];
% writematrix(Urms_limt,filename);
writematrix(combined_data,filename1);
writematrix(combined_data1,filename2);
% writematrix(PPmean_all,filename3);

combined_data2=[Umean_all1 Urms_all1];
combined_data3=[Urms_all1 Urms_power_cta];
% writematrix(Urms_limt,filename);
writematrix(combined_data2,filename3);
writematrix(combined_data3,filename4);

combined_data_all=[Umean_all Urms_all Umean_all1 Urms_all1];
writematrix(combined_data_all,filenamemax);
% writematrix(PPmean_all,filename3);
%% 

% % 応答補償
% ff=[1:NN2/2]/(NN2)*sfreq;
ff1=[1:NN2]/(NN2)*sfreq;
% ff2=ff1;
M = 0.0013189873;%input('時定数＞＞');
% 0.0004207412860 0.0004065993051  0.0004298486654 
G  = (1+2*pi*ff1*M*1i);
UU_cca = UU_CCA/Uinf;
% zz=UU_cca(:,1);
% zz_k=(fft(zz));
% ffg=zz_k.*G;


cc=0;
corrected_signal = zeros(size(UU_cca));
aafpSUMcca1 = zeros(NN2/2,1);

for kk=1:1000
    idx=[1:NN2]+(kk-1)*(sfreq);
    [idx(end) NN];
    if(idx(end)>NN)
        break;
    end
    cc=cc+1;
    % FFT → 応答補償
    X = fft(UU_cca(idx));
    Xcomp = X .* G'; % 応答補償
    % full_spectrum = [Xcomp; conj(flipud(Xcomp(2:end-1)))]; % 負の周波数も考慮 ここ
    segment = ifft(Xcomp, 'symmetric');
    
    % 補償後信号保存
    corrected_signal(idx) = segment;
    
    X2 = fft(segment)/NN2;
    Xcomp2 = X2 .* conj(X2); % 応答補償
    aafpSUMcca1 = aafpSUMcca1 +  2*Xcomp2(2:NN2/2+1);
    % % パワースペクトルの計算
    % aafpcca = Xcomp2 .* conj(Xcomp2);
    % aafpcca2_1 = 2 * aafpcca(2:NN2/2+1);
    % 
    % if kk == 1
    %     aafpSUMcca1 = aafpcca2_1;
    % else
    %     aafpSUMcca1 = aafpSUMcca1 + aafpcca2_1;
    % end
end
aafpAVEcca1 = aafpSUMcca1 / cc;

%伝達関数のチェック
figure(111)
semilogx(ff1,angle(G)/pi*180)
timtim = (0:1/sfreq*1000:40);
timtim=timtim-5;

% %応答補償後波形出力(生)(ノイズ除去前）
% figure(81)
% plot(timtim,UU_1(1:2001)/(Uinf),'k-')
% hold on
% plot(timtim,UU_CCA(1:2001)/(Uinf),'r-')
% hold on
% plot(timtim,corrected_signal(1:2001),'b--')
% xlim([10 20]);
% legend("CTA","CCA(Uncompensated)","CCA(Compensated)")
% % --- グラフの装飾 ---
% xlabel('{\it t}(s) ', 'FontSize', 15, 'FontName','Times');
% ylabel('{\it u^,/U_∞}', 'FontSize', 15, 'FontName','Times');
% legend('show', 'Location', 'best');
% box on;
% ax = gca;
% ax.XAxis.FontSize = 15;
% ax.YAxis.FontSize = 15;
% ax.XAxis.FontName = 'Times New Roman';
% ax.YAxis.FontName = ['Times New Roman'];
% ax.LineWidth = 1.3;
% set(gca,'XMinorTick','on')
% set(gca,'yMinorTick','on')
% saveas(gcf, sprintf('noisezyokyomae_hakei%ds.fig', t));

figure(4)%ノイズ除去前
loglog(ff,aafpAVE1,'k--', 'DisplayName', 'CTA');
hold on;
loglog(ff,aafpAVE,'b-', 'DisplayName', 'CCA(Uncompensated)');
loglog(ff,aafpAVEcca1, 'r--', 'DisplayName','CCA(Compensated)');
xlabel('{\it f}(Hz)')
ylabel('{\it Power}')
legend('show')
title(sprintf('CCA-CTA y=%.2f', yrmsmax));
set(gca,'Fontname','Times','Fontsize',14)
xlim([10 10^4])
ylim([10^-12 10^-2])
hold off;
% --- グラフの装飾 ---
legend('show', 'Location', 'best');
box on;
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
ax.XAxis.FontName = 'Times New Roman';
ax.YAxis.FontName = ['Times New Roman'];
ax.LineWidth = 1.3;
set(gca,'XMinorTick','on')
set(gca,'yMinorTick','on')
saveas(gcf, sprintf('noisezyokyomae_spect%ds.fig', t));

figure(937)%ノイズ除去前
loglog(ff,aafpAVE1,'k--', 'DisplayName', 'CTA', 'LineWidth',1.1);
hold on;
loglog(ff,aafpAVE,'b-', 'DisplayName', 'CCA(Uncompensated)', 'LineWidth',1.1);
loglog(ff,aafpAVEcca1, 'r--', 'DisplayName','CCA(Compensated)', 'LineWidth',1.1);
xlabel('{\it f}(Hz)')
ylabel('{\it Power}')
legend('show')
title(sprintf('CCA-CTA y=%.2f', yrmsmax));
set(gca,'Fontname','Times','Fontsize',14)
xlim([1550 2050])
ylim([10^-12 10^-2])
hold off;
% --- グラフの装飾 ---
legend('show', 'Location', 'best');
box on;
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
ax.XAxis.FontName = 'Times New Roman';
ax.YAxis.FontName = ['Times New Roman'];
ax.LineWidth = 1.3;
set(gca,'XMinorTick','on')
set(gca,'yMinorTick','on')
saveas(gcf, sprintf('noisezyokyomae_spectkakudai%ds.fig', t));
% === 高周波ノイズ除去 ===
noise = 2*10^3;%input('ノイズカットする周波数（Hz）を入力>> ');
noise_index = round(noise / (sfreq / NN)); % Hzをインデックスに変換

% FFTを取る
X_all = fft(corrected_signal);
% 正負両側をゼロ化
X_all(noise_index+1 : end-noise_index) = 0;
% 逆変換
filtered_signal = ifft(X_all, 'symmetric');

% === ノイズ除去後のFFTとスペクトル計算 ===
X_filtered = fft(filtered_signal);
aafp_filtered = X_filtered .* conj(X_filtered);
aafp_filtered2 = 2 * aafp_filtered(2:NN2/2+1);

% cc=0;
% for kk=1:1000
%     idx=[1:NN2]+(kk-1)*sfreq;%データ区間の取得
%     [idx(end) NN];
%     if(idx(end)>NN)%データ範囲が NN を超えたら終了
%         break;
%     else
%         cc=cc+1;
%         %        aa=EE(idx,1);
%         aa=filtered_signal(idx);
%         aaf=fft(aa)/NN2;%FFTして正規化NN2消していいかも%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         aafp=aaf.*conj(aaf);
%         aafp2=2*aafp(2:NN2/2+1);%パワースペクトルの対称性を考慮
%         if(kk==1)
%             aafpfilteredSUM=aafp2;
%         else
%             aafpfilteredSUM=aafpfilteredSUM+aafp2;%すべてのブロックのパワースペクトルを足し合わせていく
%         end
%     end
% end
% aafpfilteredAVE=aafpfilteredSUM/cc;%すべてのブロックのパワースペクトルを平均

% === ノイズ除去後の時間信号、スペクトルプロット ===
% %応答補償後波形出力(生)(ノイズ除去後）
% figure(815)
% plot(timtim,UU_1(1:2001)/(Uinf),'k-')
% hold on
% plot(timtim,UU_CCA(1:2001)/(Uinf),'r-')
% hold on
% plot(timtim,filtered_signal(1:2001),'b--')
% % xlim([20 30]);
% legend("CTA","CCA(Uncompensated)","CCA(Compensatedノイズ除去)")
% %pbaspect([5 1 1])
% % --- グラフの装飾 ---
% xlabel('{\it t}(s) ', 'FontSize', 15, 'FontName','Times');
% ylabel('{\it u^,/U_∞}', 'FontSize', 15, 'FontName','Times');
% legend('show', 'Location', 'best');
% box on;
% ax = gca;
% ax.XAxis.FontSize = 15;
% ax.YAxis.FontSize = 15;
% ax.XAxis.FontName = 'Times New Roman';
% ax.YAxis.FontName = ['Times New Roman'];
% ax.LineWidth = 1.3;
% set(gca,'XMinorTick','on')
% set(gca,'yMinorTick','on')
% saveas(gcf, sprintf('noisezyokyogo_hakei%d.fig', t));

figure(44)%ノイズ除去後
loglog(ff,aafpAVE1,'k--', 'DisplayName', 'CTA');
hold on;
loglog(ff,aafpAVE,'b-', 'DisplayName', 'CCA');
loglog(ff,aafp_filtered2, 'r--', 'DisplayName','CCAcompasatednoiseout');
xlabel('{\it f}(Hz)', 'FontSize', 15, 'FontName','Times');
ylabel('{\it Power}', 'FontSize', 15, 'FontName','Times');
title(sprintf('CCA-CTA y=%.2f', yrmsmax));
set(gca,'Fontname','Times','Fontsize',14)
xlim([10 10^4])
ylim([10^-12 10^-2])
hold off;
% --- グラフの装飾 ---
legend('show', 'Location', 'best');
box on;
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
ax.XAxis.FontName = 'Times New Roman';
ax.YAxis.FontName = ['Times New Roman'];
ax.LineWidth = 1.3;
set(gca,'XMinorTick','on')
set(gca,'yMinorTick','on')
saveas(gcf, sprintf('noisezyokyogo_spect%ds.fig', t));

disp('補償とノイズ除去完了');