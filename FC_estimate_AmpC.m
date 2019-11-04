function [conn_Amp, conn_Amp_Zscored, conn_Amp_static, conn_Amp_static_RegOut] = FC_estimate_AmpC(data_pre,numelec,Fs,freq,dist_electrodes)
h=waitbar(0,'Calculating Amp coupling...'); counter=0;
for i=1:numelec
    for j=1:numelec
        if i<j
            counter=counter+1;
            waitbar(counter/(numelec*(numelec-1)/2));
            conn_Amp(i,j,:)=abs(AmpCoupling(data_pre.trial{1}(i,:),data_pre.trial{1}(j,:), Fs, freq));
            conn_Amp_Zscored(i,j,:)=( conn_Amp(i,j,:)-nanmean(conn_Amp(i,j,:)) )/nanstd(conn_Amp(i,j,:));
            conn_Amp(j,i,:)=conn_Amp(i,j,:); conn_Amp_Zscored(j,i,:)=conn_Amp_Zscored(i,j,:);
        end
    end
    conn_Amp(i,i,:)=nan; conn_Amp_Zscored(i,i,:)=nan;
end
close(h)
% static FC
conn_Amp_static=nanmean(conn_Amp,3);
[conn_Amp_static_RegOut,~,~,~]=Dist_Reg_Out(conn_Amp_static,dist_electrodes);
end

