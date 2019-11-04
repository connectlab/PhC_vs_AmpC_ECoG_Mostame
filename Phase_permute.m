function [y] = Phase_permute(x)
if size(x,2)==1
    x=x';
end
signal=x(~isnan(x));
L=length(signal);
Nfft=L+(~mod(L,2));
% initializations
trialsignal_fft_phase_perm=zeros(1,Nfft);
% phase calculation
trialsignal_fft=fftshift( fft( signal,Nfft) );
trialsignal_fft_phase=angle( trialsignal_fft );
% select left half of FFT phase
ntemp=1 : floor(Nfft/2);
% permute left half
% phaseperm_temp=trialsignal_fft_phase( ntemp(randperm(length(ntemp))) ) ;
phaseperm_temp=2*pi*rand(1,length(ntemp))-pi;
% construct permuted phase
trialsignal_fft_phase_perm(ntemp(end)+2:end)=phaseperm_temp;
trialsignal_fft_phase_perm(ntemp)=-fliplr(phaseperm_temp);
% set center point phase as before
trialsignal_fft_phase_perm(ntemp(end)+1)=trialsignal_fft_phase(ntemp(end)+1);
% reconstruct permuted rawsignal
trialsignal_perm_fft=abs(trialsignal_fft).*exp(1i*trialsignal_fft_phase_perm);
trialsignal_perm=ifft( ifftshift(trialsignal_perm_fft),Nfft,'symmetric');
trialsignal_perm=trialsignal_perm(1:L);
y=x;
y(~isnan(x))=trialsignal_perm;
end

