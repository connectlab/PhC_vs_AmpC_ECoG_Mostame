%% ----------------------------------------------- Initialization
clear all
clc
% main directory
path_main='Y:\mostame2\UIUC'; cd(path_main);
% Personal functions directory
path_func='Y:\mostame2\UIUC\functions'; addpath(path_func);
% Raw data directory
path_root='Y:\mostame2\ECOGimport\'; addpath(path_root);
% Fieldtrip directory
addpath( 'Y:\mostame2\fieldtrip-20160912' );
ft_defaults;
% results path
path_results=['Y:\mostame2\UIUC' '\results\Task'];

%% ----------------------------------------------- Script
for i_sub=[1:10]
    for freq=1:5
        clc; i_sub
        freq
        Mainparham(path_main, path_root, path_results, freq,i_sub);
    end
end
%% ----------------------------------------------- Main function
function []=Mainparham(path_main, path_root, path_results, freq,i_sub)
Task='CRM'; trial_cond=1;
[subject, edata, data_pre, cfg_markeddata, electrodes_coordinate]=load_clean_data(path_root,i_sub,Task,0);
addpath(strcat('Y:\mostame2\ECOGimport\',subject)); load('BadElecs');
% Extract Data specifics
numelec=numel(edata.label);
numtrial=numel(edata.trial);
Electrodes=edata.label;
Tlim=[ edata.time{1}(1) edata.time{1}(end) ];
Fs=edata.fsample;
Time=linspace(Tlim(1), Tlim(2), size(edata.trial{1}, 2) );
%% calculate electrode distance
dist_electrodes=nan(numelec,numelec);
addpath('Y:\mostame2\ECOGimport');
fil=xlsread('MNI_Coordinates', subject);
if size(fil,2)>3
    fil(:,1:2)=[];
end
fil(BadElecs,:)=[]; fil(numelec+1:end,:)=[];
for i=1:numelec
    for j=1:numelec
        dist_electrodes(i,j)=0.1*sqrt( (fil(i,1)-fil(j,1))^2 + (fil(i,2)-fil(j,2))^2 + (fil(i,3)-fil(j,3))^2 );
    end
    dist_electrodes(i,i)=nan;
end
dist_electrodes=round(dist_electrodes,2);

%% ----------------------------------------------- Time-freq analyses configurations
Time_step=0.02;
switch freq
    case 1
        Freqrange=[5 7]; Freq_step=1; Win_length=4;
    case 2
        Freqrange=[8 13]; Freq_step=2; Win_length=6;
    case 3
        Freqrange=[14 30]; Freq_step=3; Win_length=10;
    case 4
        Freqrange=[31 60]; Freq_step=4; Win_length=20;
    case 5
        Freqrange=[61 110]; Freq_step=5; Win_length=20;
end
% set cfg1
cfg1              = [];
cfg1.output       = 'fourier';
cfg1.method       = 'mtmconvol';
cfg1.pad         = 'nextpow2';
cfg1.padtype = 'zero';
% Time-Frequency range
cfg1.toi=Tlim(1):Time_step:Tlim(end);
cfg1.foi=Freqrange(1):Freq_step:Freqrange(2);
% Time window
cfg1.t_ftimwin(:,1)=Win_length./cfg1.foi;
if freq>=3
    cfg1.taper        = 'dpss';
    % frequency smoothing
    cfg1.tapsmofrq = 0.15*cfg1.foi;
else
    cfg1.taper        = 'hanning';
end
t=cfg1.toi;
%% ----------------------------------------------- Extract power dynamics
% set cfg3
cfg3              = [];
cfg3.output       = 'pow';
cfg3.method       = 'mtmconvol';
cfg3.pad         = 'nextpow2';
cfg3.padtype = 'zero';
% Time-Frequency range
cfg3.toi=Tlim(1):Time_step:Tlim(end);
cfg3.foi=Freqrange(1):Freq_step:Freqrange(2);
% Time window
cfg3.t_ftimwin(:,1)=Win_length./cfg3.foi;
if freq>=3
    cfg3.taper        = 'dpss';
    % frequency smoothing
    cfg3.tapsmofrq = 0.15*cfg3.foi;
else
    cfg3.taper        = 'hanning';
end
% estimate short time FFT
TFR_power=ft_freqanalysis(cfg3,edata);
Power=TFR_power.powspctrm;
Power=squeeze(mean(Power,2));
% extract Z-scores with respect to pre-stimulus
for elec=1:numelec
    Power(elec,:)= ( Power(elec,:)-nanmean(Power(elec,cfg3.toi<0 & cfg3.toi>-0.5)) )/ nanstd(Power(elec,cfg3.toi<0 & cfg3.toi>-0.5));
end
%% ----------------------------------------------- Amp coupling
% estimate short time FFT
TFR=ft_freqanalysis(cfg1,edata);
% AmpAmp estimation
cfg2=[]; cfg2.removemean='yes'; cfg2.method='amplcorr'; cfg2.channelcmb= {Electrodes, Electrodes};
AmpAmp=ft_connectivityanalysis(cfg2,TFR);
FC_Amp=AmpAmp.amplcorrspctrm;
FC_Amp=squeeze(mean(FC_Amp,3)); clear temp;
FC_Amp=abs(FC_Amp);
FC_Amp_Zscored=FC_Amp;
for i=1:numelec
    for j=1:numelec
        if i<j
%             a=find( ~isnan(FC_Amp(i,j,:)),1,'first' ); b=find( ~isnan(FC_Amp(i,j,:)),1,'last' );
%             temp=squeeze(FC_Amp(i,j,:))'; temp(isnan(temp))=[];
%             [B,A]=cheby2(3,20,10/( 0.5/(cfg1.toi(2)-cfg1.toi(1)) )); temp1=filtfilt(B,A,temp); temp1=temp1/max(abs(temp1))*max(abs(temp));
%             FC_Amp(i,j,a:b)=temp;
            FC_Amp_Zscored(i,j,:)=( FC_Amp(i,j,:)-nanmean(FC_Amp(i,j,cfg1.toi<0 & cfg1.toi>-0.5)) )/nanstd(FC_Amp(i,j,cfg1.toi<0 & cfg1.toi>-0.5));
%             FC_Amp(j,i,:)=FC_Amp(i,j,:); 
            FC_Amp_Zscored(j,i,:)=FC_Amp_Zscored(i,j,:);
        end
    end
    FC_Amp(i,i,:)=nan; FC_Amp_Zscored(i,i,:)=nan;
end

%% ----------------------------------------------- PLV
% estimate PLV
cfg2=[]; cfg2.method= 'plv'; cfg2.channelcmb= {Electrodes, Electrodes};
PLV_TF=ft_connectivityanalysis(cfg2,TFR);
FC_PLV=squeeze(mean(PLV_TF.plvspctrm,2));
temp=zeros(numelec, numelec, size(FC_PLV,2));
for k=1:size(FC_PLV,1)
    i=find(strcmp(Electrodes, PLV_TF.labelcmb(k,1))>0);
    j=find(strcmp(Electrodes, PLV_TF.labelcmb(k,2))>0);
    temp(i,j,:)=FC_PLV(k,:);
end
FC_PLV=temp; clear temp; FC_PLV_Zscored=FC_PLV;
for i=1:numelec
    for j=1:numelec
        if i<j
%             a=find( ~isnan(FC_PLV(i,j,:)),1,'first' ); b=find( ~isnan(FC_PLV(i,j,:)),1,'last' );
%             temp=squeeze(FC_PLV(i,j,:))'; temp(isnan(temp))=[];
%             [B,A]=cheby2(3,20,10/( 0.5/(cfg1.toi(2)-cfg1.toi(1)) )); temp1=filtfilt(B,A,temp); temp1=temp1/max(abs(temp1))*max(abs(temp));
%             FC_PLV(i,j,a:b)=temp;
            FC_PLV_Zscored(i,j,:)=( FC_PLV(i,j,:)-nanmean(FC_PLV(i,j,cfg1.toi<0 & cfg1.toi>-0.5)) )/nanstd(FC_PLV(i,j,cfg1.toi<0 & cfg1.toi>-0.5));
%             FC_PLV(j,i,:)=FC_PLV(i,j,:); 
            FC_PLV_Zscored(j,i,:)=FC_PLV_Zscored(i,j,:);
        end
    end
    FC_PLV(i,i,:)=nan; FC_PLV_Zscored(i,i,:)=nan;
end

%% ----------------------------------------------- Exclude electrodes (Benjamini-Hochberg)
alpha=0.01;
Ignore_cond_wilcox=0; Max_cond_wilcox=0.5;
% PLV measure
t_value=[]; p_value=nan(numelec,numelec);
temp=[]; temp=randn(1,20000); temp=sort(temp,'ascend');
for i=1:numelec
    for j=1:numelec
        if i<j
            temp1=squeeze(FC_PLV_Zscored(i,j,cfg1.toi<-Ignore_cond_wilcox & cfg1.toi>-Max_cond_wilcox))';
            temp2=squeeze(FC_PLV_Zscored(i,j,cfg1.toi>Ignore_cond_wilcox & cfg1.toi<Max_cond_wilcox))';
            temp1(isnan(temp1))=[]; temp2(isnan(temp2))=[];
            s_p=nanstd([temp1 temp2]);
            t_value(i,j)=( mean(temp2)-mean(temp1) )/...
                ( s_p*(sqrt(1/length(temp1)+1/length(temp2))) );
            p_value(i,j)=1-(find(t_value(i,j)<[temp inf],1,'first')-1)/length(temp); p_value(j,i)=p_value(i,j);
        end
    end
end
temp_pval=p_value;
for i=1:numelec
    for j=1:numelec
        if i<=j
            temp_pval(i,j,:)=nan;
        end
    end
end
% right tail
clear temp; temp=temp_pval(:)'; temp=sort(temp, 'ascend'); temp(isnan(temp))=[]; 
temp=temp./[1:length(temp)]*length(temp(:));
K_hochberg=find(temp<=alpha,1,'last');
clear temp; temp=temp_pval(:)'; temp=sort(temp, 'ascend'); temp(isnan(temp))=[];
Pair_significant_PLV=double(p_value<=temp(K_hochberg));
% left tail
clear temp; temp=1-temp_pval(:)'; temp=sort(temp, 'ascend'); temp(isnan(temp))=[];
temp=temp./[1:length(temp)]*length(temp(:));
K_hochberg=find(temp<=alpha,1,'last');
clear temp; temp=1-temp_pval(:)'; temp=sort(temp, 'ascend'); temp(isnan(temp))=[];
Pair_significant_PLV=Pair_significant_PLV - double( (1-p_value)<=temp(K_hochberg) );

% Amp measure
t_value=[]; p_value=nan(numelec,numelec);
temp=[]; temp=randn(1,20000); temp=sort(temp,'ascend');
for i=1:numelec
    for j=1:numelec
        if i<j
            temp1=squeeze(FC_Amp_Zscored(i,j,cfg1.toi<-Ignore_cond_wilcox & cfg1.toi>-Max_cond_wilcox))';
            temp2=squeeze(FC_Amp_Zscored(i,j,cfg1.toi>Ignore_cond_wilcox & cfg1.toi<Max_cond_wilcox))';
            temp1(isnan(temp1))=[]; temp2(isnan(temp2))=[];
            s_p=nanstd([temp1 temp2]);
            t_value(i,j)=( mean(temp2)-mean(temp1) )/...
                ( s_p*(sqrt(1/length(temp1)+1/length(temp2))) );
            p_value(i,j)=1-(find(t_value(i,j)<[temp inf],1,'first')-1)/length(temp); p_value(j,i)=p_value(i,j);
        end
    end
end
temp_pval=p_value;
for i=1:numelec
    for j=1:numelec
        if i<=j
            temp_pval(i,j,:)=nan;
        end
    end
end
% right tail
clear temp; temp=temp_pval(:)'; temp=sort(temp, 'ascend'); temp(isnan(temp))=[]; 
temp=temp./[1:length(temp)]*length(temp(:));
K_hochberg=find(temp<=alpha,1,'last');
clear temp; temp=temp_pval(:)'; temp=sort(temp, 'ascend'); temp(isnan(temp))=[];
Pair_significant_Amp=double(p_value<=temp(K_hochberg));
% left tail
clear temp; temp=1-temp_pval(:)'; temp=sort(temp, 'ascend'); temp(isnan(temp))=[];
temp=temp./[1:length(temp)]*length(temp(:));
K_hochberg=find(temp<=alpha,1,'last');
clear temp; temp=1-temp_pval(:)'; temp=sort(temp, 'ascend'); temp(isnan(temp))=[];
Pair_significant_Amp=Pair_significant_Amp - double( (1-p_value)<=temp(K_hochberg) );

% merge two measures
Pair_significant=Pair_significant_PLV | Pair_significant_Amp;
for i=1:numelec
    Pair_significant_PLV(i,i)=0; Pair_significant_Amp(i,i)=0; Pair_significant(i,i)=0;
end

%% ----------------------------------------------- Select trial period of interest
switch trial_cond % if 0, check what matrices are you using: Zscored or no!
    case 0 % Pre-stimulus
        conn_PLV=FC_PLV(:,:,cfg1.toi<0); conn_Amp=FC_Amp(:,:,cfg1.toi<0);
        Toi=cfg1.toi(cfg1.toi<0);
        Tlim_cond=[Tlim(1) 0];
    case 1 % Post-stimulus
        conn_PLV=FC_PLV_Zscored(:,:,cfg1.toi>0); conn_Amp=FC_Amp_Zscored(:,:,cfg1.toi>0);
        Toi=cfg1.toi(cfg1.toi>0);
        Tlim_cond=[0 Tlim(2)];
    case 2 % Whole trial
        conn_PLV=FC_PLV_Zscored; conn_Amp=FC_Amp_Zscored;
        Toi=cfg1.toi;
        Tlim_cond=Tlim;
end
%% ----------------------------------------------- Generate null PDF of maximum correlations
Max_Lag=0.5; %ms
Max_Lag=Max_Lag/(Toi(2)-Toi(1));
R=500;
Max_corr_surr=nan(numelec,numelec,R); Zerolag_corr_surr=nan(numelec,numelec,R);
lag_diff_nearest_surr=nan(numelec,numelec);
Threshold_corr_pos=nan(numelec,numelec); Threshold_corr_neg=nan(numelec,numelec);
h=waitbar(0, 'Please wait...Generating surrogate correlation values...'); counter=0;
for i=1:numelec
    for j=1:numelec
        if i<j && Pair_significant(i,j)
            counter=counter+1;
            waitbar( counter/ (0.5*sum(sum(Pair_significant~=0))) );
            L=length(find(t>0));
            for repeat=1:R
                temp1=squeeze(FC_PLV_Zscored(i,j,:))';
                temp1(isnan(temp1))=[];
                temp1=Phase_permute(temp1);
                temp1=temp1(end-L+1:end);
                temp2=squeeze(FC_Amp_Zscored(i,j,:))';
                temp2(isnan(temp2))=[];
                temp2=Phase_permute(temp2);
                temp2=temp2(end-L+1:end);
                %calculate correlation
                [acorr, lag]=xcorr(temp1-nanmean(temp1), temp2-nanmean(temp2),'coeff');
                Zerolag_corr_surr(i,j,repeat)=acorr(lag==0); Zerolag_corr_surr(j,i,repeat)=Zerolag_corr_surr(i,j,repeat);
%                 % find out lag
%                 a=round(Max_Lag); b=1+floor(length(lag)/2);
%                 lag=lag(b-a:b+a); acorr=acorr(b-a:b+a);
%                 dt_acorr=[diff(acorr) 0];
%                 temp=[]; temp=dt_acorr.*circshift(dt_acorr,1); tmp=find(temp<0)-ceil(0.5*length(lag));
%                 lag_diff_nearest_surr(i,j)=tmp( find(abs(tmp)==min(abs(tmp)),1,'first') );
%                 lag_diff_nearest_surr(j,i)=lag_diff_nearest_surr(i,j);
%                 % find out corresponding peak
%                 Max_corr_surr(i,j,repeat)=acorr(lag==lag_diff_nearest_surr(i,j));
%                 Max_corr_surr(j,i,repeat)=Max_corr_surr(i,j,repeat);
            end
        end
    end
end
close(h)

%% ----------------------------------------------- estimate maximum correlations and lags
% estimate correlations
Max_corr=nan(numelec,numelec); Zerolag_corr=nan(numelec,numelec);
lag_diff_nearest=nan(numelec,numelec);
for i=1:numelec
    for j=1:numelec
        if i<j && Pair_significant(i,j)
            temp1=[]; temp2=[]; temp1=squeeze(conn_PLV(i,j,:))'; temp2=squeeze(conn_Amp(i,j,:))';
            temp1(isnan(temp1))=[]; temp2(isnan(temp2))=[];
            [acorr, lag]=xcorr(temp1-nanmean(temp1), temp2-nanmean(temp2),'coeff');
            Zerolag_corr(i,j)=acorr(lag==0); Zerolag_corr(j,i)=Zerolag_corr(i,j);
%             % find out best lag
%             a=round(Max_Lag); b=1+floor(length(lag)/2);
%             lag=lag(b-a:b+a); acorr=acorr(b-a:b+a);
%             dt_acorr=[diff(acorr) 0];
%             temp=[]; temp=dt_acorr.*circshift(dt_acorr,1); tmp=find(temp<0)-ceil(0.5*length(lag));
%             lag_diff_nearest(i,j)=tmp( find(abs(tmp)==min(abs(tmp)),1,'first') );
%             lag_diff_nearest(j,i)=lag_diff_nearest(i,j);
%             % find out corresponding peak
%             Max_corr(i,j)=acorr(lag==lag_diff_nearest(i,j)); Max_corr(j,i)=Max_corr(i,j);
        end
    end
end
% --------------------------- test
Max_corr=Zerolag_corr;
% correct lag matrix
% lag_diff_nearest=lag_diff_nearest.*(Toi(2)-Toi(1))*Fs; lag_diff_nearest(Pair_significant~=1)=nan;
% filter large lags
% Pair_significant(abs(lag_diff_nearest)>250)=0;
% lag_diff_nearest(Pair_significant~=1)=nan;
Max_corr(Pair_significant~=1)=nan; Zerolag_corr(Pair_significant~=1)=nan;
%% fisher transformation
% for i=1:size(Max_corr,1)
%     for j=1:size(Max_corr,1)
%         Max_corr_fisher(i,j)=0.5*log( (1+Max_corr(i,j))/(1-Max_corr(i,j)) );
%         Max_corr_fisher(j,i)=Max_corr_fisher(i,j);
%     end
% end
% for k=1:size(Zerolag_corr_surr,3)
%     for i=1:size(Zerolag_corr_surr,1)
%         for j=1:size(Zerolag_corr_surr,1)
%             Max_corr_surr_fisher(i,j,k)=0.5*log( (1+Zerolag_corr_surr(i,j,k))/(1-Zerolag_corr_surr(i,j,k)) );
%             Max_corr_surr_fisher(j,i,k)=Max_corr_surr_fisher(i,j,k);
%         end
%     end
% end
%% ----------------------------------------------- FDR control statistical test
% MCP Hochberg FDR for positive values
alpha=0.05;
pval=nan(numelec,numelec);
for i=1:numelec
    for j=1:numelec
        if i<j && Pair_significant(i,j)==1
            clear temp; temp=sort( squeeze(Zerolag_corr_surr(i,j,:)), 'ascend' )';
            pval(i,j)=1- (find(Max_corr(i,j)<[temp inf],1,'first')-1)/length(temp); pval(j,i)=pval(i,j);
        end
    end
end
temp_pval=pval;
for i=1:numelec
    for j=1:numelec
        if i<=j
            temp_pval(i,j)=nan;
        end
    end
end
% positive tail
clear temp; temp=sort(temp_pval(:), 'ascend')'; temp(isnan(temp))=[];
temp=temp./[1:length(temp)]*length(temp);
K_hochberg=find(temp<=alpha,1,'last');
if ~isempty(K_hochberg)
    clear temp; temp=sort(temp_pval(:), 'ascend')'; temp(isnan(temp))=[];
    Max_corr_significant_hochberg=double(pval<=temp(K_hochberg));
else
    Max_corr_significant_hochberg=zeros(size(pval));
end
% negative tail
clear temp; temp=sort(1-temp_pval(:), 'ascend')'; temp(isnan(temp))=[];
temp=temp./[1:length(temp)]*length(temp);
K_hochberg=find(temp<=alpha,1,'last');
if ~isempty(K_hochberg)
    clear temp; temp=sort(1-temp_pval(:), 'ascend')'; temp(isnan(temp))=[];
    Max_corr_significant_hochberg=Max_corr_significant_hochberg -double((1-pval)<=temp(K_hochberg));
end
%exclude irrelevant electrodes
Max_corr_significant_hochberg(Pair_significant==0)=nan;

%% --------------------------------------- estimate FC dynamics from continuous data
%% clean continuous data (correct time jumps on this data)
for i=size(cfg_markeddata.artfctdef.visual.artifact,1):-1:1
    data_pre.trial{1}(:, cfg_markeddata.artfctdef.visual.artifact(i,1):cfg_markeddata.artfctdef.visual.artifact(i,2) )=[];
end
data_pre.time{1}=(0:length(data_pre.trial{1}(1,:))-1 )./Fs;
data_pre.sampleinfo=length( data_pre.time{1} );

%% estimate AmpC for continuous data
conn_Amp_cont=[]; conn_Amp_Zscored_cont=[];
[conn_Amp_cont, conn_Amp_Zscored_cont, conn_Amp_static, conn_Amp_static_RegOut]=FC_estimate_AmpC(data_pre,numelec,Fs,freq,dist_electrodes);
%% estimate PhC for continuous data
conn_PLV_cont=[]; conn_PLV_Zscored_cont=[];
[conn_PLV_cont, conn_PLV_Zscored_cont, conn_PLV_static, conn_PLV_static_RegOut]=FC_estimate_PLV(data_pre,numelec,Fs,freq,dist_electrodes);
%% find out length of time window
% for PhC
A=[]; B=[];
for i=1:numelec
    for j=1:numelec
        if i<j
            temp=[];
            temp=squeeze(conn_PLV_Zscored_cont(i,j,:))'; temp=temp(~isnan(temp));
            [a w]=cpsd(temp,temp,[],[],1024); w=w/(2*pi);
            % find right side FWHM
            temp=[]; temp=abs((a-0.5*max(a))); temp(w<w(a==max(a)))=nan; temp=w( temp-nanmin(temp) == min( temp-nanmin(temp) ) );
%             plot(w,a,'linewidth',2); line([temp temp],[0 a(w==temp)+2], 'color','k')
            A=[A w(a==max(a))];
            B=[B temp];
        end
    end
end
W1=median(A); W1_FWHM=median(B);
% for AmpC
A=[]; B=[];
for i=1:numelec
    for j=1:numelec
        if i<j
            temp=[];
            temp=squeeze(conn_Amp_Zscored_cont(i,j,:))'; temp=temp(~isnan(temp));
            [a w]=cpsd(temp,temp,[],[],1024); w=w/(2*pi);
            % find right side FWHM
            temp=[]; temp=abs((a-0.5*max(a))); temp(w<w(a==max(a)))=nan; temp=w( temp-nanmin(temp) == min( temp-nanmin(temp) ) );
%             plot(w,a,'linewidth',2); line([temp temp],[0 a(w==temp)+2], 'color','k')
            A=[A w(a==max(a))];
            B=[B temp];
        end
    end
end
W2=median(A); W2_FWHM=median(B);
Corr_step=round(1/min(W1, W2),1,'significant');
if Corr_step>100
    fprintf('<<<<<<<<<<<<<<<WARNING! Corr_Step is longer than 100s = %d>>>>>>>>>>>>>>\n',Corr_step)
    Corr_step=100;
elseif Corr_step<20
    fprintf('<<<<<<<<<<<<<<<WARNING! Corr_Step is shorter than 20s = %d>>>>>>>>>>>>>>\n',Corr_step)
    Corr_step=20;
end
clear A B a w W1 W2
%% assess correlation between the two dynamics
h=waitbar(0,'Calculating dynamic correlations...'); counter=0;
Corr=nan(size(conn_PLV_static));
Corr_alltime=[];
for k=1:Corr_step/2:(size(conn_Amp_cont,3)-Corr_step+1)
    counter=counter+1;
    waitbar(counter/(numel(1:Corr_step/2:size(conn_Amp_cont,3)-Corr_step+1)));
    for i=1:numelec
        for j=1:numelec
            if i<j && Pair_significant(i,j)==1
                temp1=[]; temp1=conn_PLV_Zscored_cont(i,j,k:k+Corr_step-1); temp1=squeeze(temp1)'; temp1(isnan(temp1))=[];
                temp2=[]; temp2=conn_Amp_Zscored_cont(i,j,k:k+Corr_step-1); temp2=squeeze(temp2)'; temp2(isnan(temp2))=[];
                [acorr lag]=xcorr(temp1-nanmean(temp1),temp2-nanmean(temp2),'coeff');
                Corr(i,j)=acorr(lag==0); Corr(j,i)=Corr(i,j);
            end
        end
        Corr(i,i)=nan;
    end
    Corr_alltime=cat(3,Corr_alltime,Corr);
end
close(h); clear h;
%% permute FC dynamics for statistical test
R=500;
Corr_surr_repeat=nan([size(Corr_alltime), R]);
h=waitbar(0,'generating  correlations...'); counter=0;
for repeat=1:R
    counter=counter+1;
    waitbar(counter/R);
    Corr_surr_alltime=nan(size(Corr_alltime));
    for i=1:numelec
        for j=1:numelec
            if i<j && Pair_significant(i,j)==1
                temp=[]; temp=squeeze(conn_PLV_Zscored_cont(i,j,:))'; temp1=Phase_permute(temp);
                temp=[]; temp=squeeze(conn_Amp_Zscored_cont(i,j,:))'; temp2=Phase_permute(temp);
                count=0;
                for k=1:Corr_step/2:(size(conn_Amp_cont,3)-Corr_step+1)
                    count=count+1;
                    tmp1=temp1(k:k+Corr_step-1); tmp2=temp2(k:k+Corr_step-1); 
                    tmp1(isnan(tmp1))=[]; tmp2(isnan(tmp2))=[];
                    [acorr lag]=xcorr(tmp1-nanmean(tmp1),tmp2-nanmean(tmp2),'coeff');
                    Corr_surr_alltime(i,j,count)=acorr(lag==0); Corr_surr_alltime(j,i,count)=Corr_surr_alltime(i,j,count);
                end
            end
        end
    end
    Corr_surr_repeat(:,:,:,repeat)=Corr_surr_alltime;
end
close(h); clear h PLV_perm Amp_perm

%% ----------------------------------------------- FDR control statistical test
% MCP Hochberg FDR for positive values
alpha=0.05;
pval=nan(size(Corr_alltime));
for k=1:size(Corr_alltime,3)
    for i=1:numelec
        for j=1:numelec
            if i<j && Pair_significant(i,j)==1
                clear temp; temp=sort( squeeze(Corr_surr_repeat(i,j,k,:)), 'ascend' )';
                pval(i,j,k)=1- (find(Corr_alltime(i,j,k)<[temp inf],1,'first')-1)/length(temp); pval(j,i,k)=pval(i,j,k);
            end
        end
    end
end
temp_pval=pval;
for i=1:numelec
    for j=1:numelec
        if i<=j
            temp_pval(i,j,:)=nan;
        end
    end
end
% positive tail
clear temp; temp=sort(temp_pval(:), 'ascend')'; temp(isnan(temp))=[];
temp=temp./[1:length(temp)]*length(temp);
K_hochberg=find(temp<=alpha,1,'last');
if ~isempty(K_hochberg)
    clear temp; temp=sort(temp_pval(:), 'ascend')'; temp(isnan(temp))=[];
    Corr_alltime_significant_hochberg=double(pval<=temp(K_hochberg));
else
    Corr_alltime_significant_hochberg=zeros(size(pval));
end
% negative tail
clear temp; temp=sort(1-temp_pval(:), 'ascend')'; temp(isnan(temp))=[];
temp=temp./[1:length(temp)]*length(temp);
K_hochberg=find(temp<=alpha,1,'last');
if ~isempty(K_hochberg)
    clear temp; temp=sort(1-temp_pval(:), 'ascend')'; temp(isnan(temp))=[];
    Corr_alltime_significant_hochberg=Corr_alltime_significant_hochberg -double((1-pval)<=temp(K_hochberg));
end
%exclude irrelevant electrodes
for i=1:size(Corr_alltime_significant_hochberg,3)
    temp=squeeze(Corr_alltime_significant_hochberg(:,:,i)); temp(Pair_significant==0)=nan;
    Corr_alltime_significant_hochberg(:,:,i)=temp;
end

%% ----------------------------------------------- extract correlation degrees
corr_degree_pos=nanmean(Max_corr.*(Max_corr_significant_hochberg>0),2);
corr_degree_neg=nanmean(Max_corr.*(Max_corr_significant_hochberg<0),2); corr_degree_neg=abs(corr_degree_neg);
%% ----------------------------------------------- static connectivity
% ------------ Baseline
cfg=[]; cfg.toilim=[Tlim(1) 0];
edata_trim=ft_redefinetrial(cfg, edata);
% set cfg1
cfg0              = [];
cfg0.output       = 'fourier';
cfg0.method       = 'mtmfft';
cfg0.pad         = 'nextpow2';
cfg0.padtype = 'zero';
cfg0.foi=Freqrange(1):Freq_step:Freqrange(2);
if freq>=3
    cfg0.taper        = 'dpss';
    % frequency smoothing
    cfg0.tapsmofrq = 5;
else
    cfg0.taper        = 'hanning';
end
% estimate short time FFT
TFR=ft_freqanalysis(cfg0,edata_trim);
% AmpAmp estimation
cfg2=[]; cfg2.removemean='yes'; cfg2.method='amplcorr'; cfg2.channelcmb= {Electrodes, Electrodes};
AmpAmp=ft_connectivityanalysis(cfg2,TFR);
FC_Amp_static_baseline=AmpAmp.amplcorrspctrm;
FC_Amp_static_baseline=squeeze(mean(FC_Amp_static_baseline,3)); clear temp;
FC_Amp_static_baseline=abs(FC_Amp_static_baseline);
for elec=1:numelec
    FC_Amp_static_baseline(elec,elec)=nan;
end
% PLV estimation
cfg2=[]; cfg2.method= 'plv'; cfg2.channelcmb= {Electrodes, Electrodes};
PLV_TF=ft_connectivityanalysis(cfg2,TFR);
FC_PLV_static_baseline=squeeze(mean(PLV_TF.plvspctrm,2));
temp=zeros(numelec, numelec);
for k=1:size(FC_PLV_static_baseline,1)
    i=find(strcmp(Electrodes, PLV_TF.labelcmb(k,1))>0);
    j=find(strcmp(Electrodes, PLV_TF.labelcmb(k,2))>0);
    temp(i,j)=FC_PLV_static_baseline(k);
end
FC_PLV_static_baseline=temp; clear temp;
for elec=1:numelec
    FC_PLV_static_baseline(elec,elec)=nan;
end
%% regress out distance
[FC_PLV_static_RegOut_baseline,~,~,~]=Dist_Reg_Out(FC_PLV_static_baseline,dist_electrodes);
[FC_Amp_static_RegOut_baseline,~,~,~]=Dist_Reg_Out(FC_Amp_static_baseline,dist_electrodes);
% time-averaged static FC
FC_PLV_static_timeaveraged_baseline=nanmean(FC_PLV(:,:, t<0),3);
FC_Amp_static_timeaveraged_baseline=nanmean(FC_Amp(:,:, t<0),3);
% regress out distance
[FC_PLV_static_timeaveraged_RegOut_baseline,~,~,~]=Dist_Reg_Out(FC_PLV_static_timeaveraged_baseline,dist_electrodes);
[FC_Amp_static_timeaveraged_RegOut_baseline,~,~,~]=Dist_Reg_Out(FC_Amp_static_timeaveraged_baseline,dist_electrodes);

% ------------ post-stimlus
cfg=[]; cfg.toilim=[0 Tlim(2)];
edata_trim=ft_redefinetrial(cfg, edata);
% set cfg1
cfg0              = [];
cfg0.output       = 'fourier';
cfg0.method       = 'mtmfft';
cfg0.pad         = 'nextpow2';
cfg0.padtype = 'zero';
cfg0.foi=Freqrange(1):Freq_step:Freqrange(2);
if freq>=3
    cfg0.taper        = 'dpss';
    % frequency smoothing
    cfg0.tapsmofrq = 5;
else
    cfg0.taper        = 'hanning';
end
% estimate short time FFT
TFR=ft_freqanalysis(cfg0,edata_trim);
% AmpAmp estimation
cfg2=[]; cfg2.removemean='yes'; cfg2.method='amplcorr'; cfg2.channelcmb= {Electrodes, Electrodes};
AmpAmp=ft_connectivityanalysis(cfg2,TFR);
FC_Amp_static=AmpAmp.amplcorrspctrm;
FC_Amp_static=squeeze(mean(FC_Amp_static,3)); clear temp;
FC_Amp_static=abs(FC_Amp_static);
for elec=1:numelec
    FC_Amp_static(elec,elec)=nan;
end
% PLV estimation
cfg2=[]; cfg2.method= 'plv'; cfg2.channelcmb= {Electrodes, Electrodes};
PLV_TF=ft_connectivityanalysis(cfg2,TFR);
FC_PLV_static=squeeze(mean(PLV_TF.plvspctrm,2));
temp=zeros(numelec, numelec);
for k=1:size(FC_PLV_static,1)
    i=find(strcmp(Electrodes, PLV_TF.labelcmb(k,1))>0);
    j=find(strcmp(Electrodes, PLV_TF.labelcmb(k,2))>0);
    temp(i,j)=FC_PLV_static(k);
end
FC_PLV_static=temp; clear temp;
for elec=1:numelec
    FC_PLV_static(elec,elec)=nan;
end
% regress out distance
[FC_PLV_static_RegOut,~,~,~]=Dist_Reg_Out(FC_PLV_static,dist_electrodes);
[FC_Amp_static_RegOut,~,~,~]=Dist_Reg_Out(FC_Amp_static,dist_electrodes);
% time-averaged static FC
FC_PLV_static_timeaveraged=nanmean(FC_PLV(:,:, t>0),3);
FC_Amp_static_timeaveraged=nanmean(FC_Amp(:,:, t>0),3);
% regress out distance
[FC_PLV_static_timeaveraged_RegOut,~,~,~]=Dist_Reg_Out(FC_PLV_static_timeaveraged,dist_electrodes);
[FC_Amp_static_timeaveraged_RegOut,~,~,~]=Dist_Reg_Out(FC_Amp_static_timeaveraged,dist_electrodes);
% surrogate static FC matrices
for repeat=1:1000
    FC_PLV_static_regout_surr(:,:,repeat)=Phase_permute_2D(FC_PLV_static_RegOut);
    FC_Amp_static_regout_surr(:,:,repeat)=Phase_permute_2D(FC_Amp_static_RegOut);
    FC_PLV_static_surr(:,:,repeat)=Phase_permute_2D(FC_PLV_static);
    FC_Amp_static_surr(:,:,repeat)=Phase_permute_2D(FC_Amp_static);
end
clear repeat
%% ----------------------------------------------- save data
clear temp
cd(path_results);
txt=sprintf('AmpPhCorr_%s_Cond%d_Freq%dto%d_smooth_abs.mat',subject,trial_cond,Freqrange(1),Freqrange(2));
save(txt)
cd(path_main);
end
