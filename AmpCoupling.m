function [Amp] = AmpCoupling(x,y,Fs, freq)
% reduce computational cost
if Fs>250
    x=downsample(x,floor(Fs/250)); y=downsample(y,floor(Fs/250));
    Fs=Fs/floor(Fs/250);
end
% trim extra data
L=min(length(x),length(y));
if length(x)>L
    x(L+1:end)=[];
elseif length(y)>L
    y(L+1:end)=[];
end
% initialize time vector
t=0:1/Fs:(L-1)/Fs;
Time_step=1;
switch freq
    case 1
        Freqrange=[5 7]; Win_length=75;
    case 2
        Freqrange=[8 13]; Win_length=100;
    case 3
        Freqrange=[14 30]; Win_length=200;
    case 4
        Freqrange=[31 60]; Win_length=400;
    case 5
        Freqrange=[61 110]; Win_length=800;
end
Win_step=floor(Win_length/mean(Freqrange));
% band-pass filter
[B,A]=cheby2(4,40,Freqrange/Fs*2, 'bandpass');
x_bp=filtfilt(B,A,x); x_bp=x_bp/max(abs(x_bp))*max(abs(x));
y_bp=filtfilt(B,A,y); y_bp=y_bp/max(abs(y_bp))*max(abs(y));
% estimate FC in each shifting time window
starttime=Win_step/2;
while starttime + Win_step/2 <=t(end)
    x_bp_trim=x_bp( find(t>=starttime-Win_step/2,1,'first') : find(t<starttime+Win_step/2,1,'last'));
    y_bp_trim=y_bp( find(t>=starttime-Win_step/2,1,'first') : find(t<starttime+Win_step/2,1,'last'));
    x_bp_trim_env=abs(hilbert(x_bp_trim));
    y_bp_trim_env=abs(hilbert(y_bp_trim));
    [acorr, lag]=xcorr(x_bp_trim_env-mean(x_bp_trim_env), y_bp_trim_env-mean(y_bp_trim_env), 'coeff');
    Amp( 1+ (starttime-Win_step/2)/Time_step )=acorr(lag==0);
    starttime=starttime+Time_step;
end
% smooth dynamic
% [B,A]=cheby2(3,20,2*0.3*Time_step); temp1=filtfilt(B,A,Amp); temp1=temp1/max(abs(temp1))*max(abs(Amp));
% Amp=temp1; clear temp1;
% add nan values to the tails
temp=[]; temp=floor(t(end))/Time_step;
temp=temp-length(Amp);
Amp=[nan(1,floor(temp/2)) Amp nan(1,ceil(temp/2))];
end
