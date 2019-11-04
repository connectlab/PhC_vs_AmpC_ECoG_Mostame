function [FC_PLV_RegOut,Xq,Vq,Vq_RegOut] = Dist_Reg_Out(FC_PLV,dist_electrodes)
FC_PLV=FC_PLV(:)'; dist_electrodes=dist_electrodes(:)';
% find out samples
x=0:1.5:15;
Samples=zeros(1,numel(x));
for i=1:numel(x)
    temp=FC_PLV( dist_electrodes>x(i)-(x(2)-x(1))/2 & dist_electrodes<x(i)+(x(2)-x(1))/2 );
    temp(isnan(temp))=[];
    Samples(i)=nanmean(temp(:));
end
Xq=x(1):0.01:x(end); Xq=round(Xq,2);
Vq=interp1(x,Samples,Xq,'PCHIP');
Vq(Xq>max(dist_electrodes))=[]; Xq(Xq>max(dist_electrodes))=[];
Vq(Xq<min(dist_electrodes))=[]; Xq(Xq<min(dist_electrodes))=[];
% subtract
dist_electrodes=round(dist_electrodes,2);
for i=1:numel(dist_electrodes)
    tmp=find(dist_electrodes(i)==Xq,2);
    if ~isnan(tmp)
        FC_PLV_RegOut( i )=FC_PLV( i ) -Vq( tmp );
    else
        FC_PLV_RegOut( i )=nan;
    end
end
FC_PLV_RegOut=reshape(FC_PLV_RegOut,sqrt(numel(FC_PLV_RegOut)), sqrt(numel(FC_PLV_RegOut)));
Samples=zeros(1,numel(x));
for i=1:numel(x)
    temp=FC_PLV_RegOut(dist_electrodes>x(i)-1 & dist_electrodes<x(i)+1);
    temp(isnan(temp))=[];
    Samples(i)=nanmean(temp(:));
end
Vq_RegOut=interp1(x,Samples,Xq,'PCHIP');
Vq_RegOut(Xq>max(dist_electrodes))=[]; Xq(Xq>max(dist_electrodes))=[];
Vq_RegOut(Xq<min(dist_electrodes))=[]; Xq(Xq<min(dist_electrodes))=[];
end

