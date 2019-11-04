function [y] = Phase_permute_2D(x)
A=isnan(x);
x(A)=0;
x_fft=fftshift(fft2(x));
x_fft_angle=angle(x_fft);
temp=x_fft_angle(1:1+floor(0.5*size(x_fft_angle,1)),:); [s1 s2]=size(temp);
temp=temp(randperm(numel(temp))); temp=reshape(temp,s1,s2);
x_fft_angle_perm=x_fft_angle;
x_fft_angle_perm(1:s1,:)=temp;
temp=-fliplr(flipud(temp)); temp(1:2*s1-size(x,1),:)=[];
x_fft_angle_perm(1+s1:end,:)=temp;
x_fft_angle_perm(s1,s1)=0;
y_fft=abs(x_fft).*exp(1i*x_fft_angle_perm);
y=ifft2(fftshift(y_fft),'symmetric');
y(A)=nan;
end


