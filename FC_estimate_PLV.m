function [conn_PLV, conn_PLV_Zscored, conn_PLV_static, conn_PLV_static_RegOut] = FC_estimate_PLV(data_pre,numelec,Fs,freq,dist_electrodes)
h=waitbar(0,'Calculating phase coupling (PLV)...'); counter=0;
for i=1:numelec
    for j=1:numelec
        if i<j
            counter=counter+1;
            waitbar(counter/(numelec*(numelec-1)/2));
            conn_PLV(i,j,:)=PLV_Sepideh(data_pre.trial{1}(i,:),data_pre.trial{1}(j,:), Fs, freq);
            conn_PLV_Zscored(i,j,:)=( conn_PLV(i,j,:)-nanmean(conn_PLV(i,j,:)) )/nanstd(conn_PLV(i,j,:));
            conn_PLV(j,i,:)=conn_PLV(i,j,:); conn_PLV_Zscored(j,i,:)=conn_PLV_Zscored(i,j,:);
        end
    end
    conn_PLV(i,i,:)=nan; conn_PLV_Zscored(i,i,:)=nan;
end
close(h)
% static FC
conn_PLV_static=nanmean(conn_PLV,3);
[conn_PLV_static_RegOut,~,~,~]=Dist_Reg_Out(conn_PLV_static,dist_electrodes);
end

