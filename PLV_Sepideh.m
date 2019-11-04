function [PLV] = PLV_Sepideh(x,y,Fs, freq)
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
t=0:1/Fs:L/Fs;
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
while starttime + Win_step/2<=t(end)
    x_bp_trim=x_bp( find(t>starttime -Win_step/2,1,'first') : find(t>starttime+Win_step/2,1,'first') );
    y_bp_trim=y_bp( find(t>starttime -Win_step/2,1,'first') : find(t>starttime+Win_step/2,1,'first') );
    % extract raw phases
    Phases=[];
    Phases(1,:)=angle(hilbert(x_bp_trim)); Phases(2,:)=angle(hilbert(y_bp_trim));
    % smoothen phases
    for elec=1:2
        ph=Phases(elec,:);
        for time=1:length(ph)-1
            if ph(time+1)-ph(time)>2.4
                ph(time+1:end)=ph(time+1:end)-2*pi;
            elseif ph(time+1)-ph(time)<-2.4
                ph(time+1:end)=ph(time+1:end)+2*pi;
            end
        end
        Phases(elec,:)=ph;
    end
    clear ph
    PLV( 1+ (starttime-Win_step/2)/Time_step )=abs(mean(exp( 1i*(Phases(1,:)-Phases(2,:)) )));
    starttime=starttime+Time_step;
end
% smooth curves
% [B,A]=cheby2(3,20,2*0.3*Time_step); temp1=filtfilt(B,A,PLV); temp1=temp1/max(abs(temp1))*max(abs(PLV));
% PLV=temp1; clear temp1;
% add nan values to the tails
temp=[]; temp=floor(t(end))/Time_step;
temp=temp-length(PLV);
PLV=[nan(1,floor(temp/2)) PLV nan(1,ceil(temp/2))];
end
