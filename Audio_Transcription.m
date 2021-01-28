%% Part II: transcript the music notes from some audio input recording

clc;
close all;
clear all;
[x1,fs1]=audioread('ech24.wav');
x=resample(x1,8000,44100);   %to sample at fs=8000Hz
fs1
fs=8000;
Ts=1/fs;
f1=figure;
t=0:1/fs:(length(x)-1)/fs;   %time-axis
plot(t,x);                   %plot in time domain
xlabel('Time');
ylabel('Amplitude');
title('original signal-time domain');

f=0:fs/(length(x)-1):fs;     %frequency axis
xfft=abs(fft(x));            %taking FFT of the signal
f2=figure;
plot(f,xfft);           %plot in frequency domain
xlabel('Frequency');
ylabel('FFT');
title('Frequency domain representation'); 

%% Low Pass Filtering  
%using Butterworth LPF

Fp=500;               % passband frequency
Fs=1000;              % stopband frequency
Wp=2*pi*Fp/fs;
Ws=2*pi*Fs/fs;
Wp1=2/Ts*tan(Wp/2);
Ws1=2/Ts*tan(Ws/2);
rp=1;                 % passband ripple
rs=30;                % stopband attenuation
[n,Wn]=buttord(Wp1,Ws1,rp,rs,'s');   % finds the minimum order n and cutoff frequencies Wn
[Z,P,K]=buttap(n);     %butterworth filter prototype
[B,A]=zp2tf(Z,P,K);    %to get polynomial transfer function representation
[b,a]=lp2lp(B,A,Wn);   % to get a lowpass filter with cutoff angular frequency
[bz,az]=bilinear(b,a,fs);
[H,W]=freqz(bz,az);    %frequency response of the filter
y = filter(bz,az,x);   %to filter x with coefficients a and b
f3=figure;
plot(y); 
xlabel('Time');
ylabel('Amplitude');
title('Filtered signal - time domain');%time domain
yfft=abs(fft(y)); 
f4=figure;
plot(f,yfft);         %frequency domain
xlabel('Frequency');
ylabel('FFT');
title('Filtered signal- frequency domain');

filename='ech24_filtered.wav';
audiowrite(filename,y,fs1);

%% time-frequency analysis

N =512;         %length of the window
h=hamming(N);
Freq_scale=0:1:1000;     %frequency scale
[S,Freq_scale,t1]=spectrogram(y(:,1),h,N-12,Freq_scale,fs,'yaxis');
f5=figure;
S_fft=abs(S);
imagesc(t1,Freq_scale,S_fft);   %to display the data as an image
xlabel('Time');
ylabel('Frequency');

%% Note transcription

[val,ind]=max(S_fft); %to get the index of max value
maxi = max(ind);
mini = min(ind);

i=1;
for j=1:2534                    
    %the range is specified with reference to the frequency table(octave 3)
    if ind(i)>255 && ind(i)<265
        n(i)=1;
    elseif ind(i)>270 && ind(i)<285
        n(i)=2;
    elseif ind(i)>286 && ind(i)<300
        n(i)=3;   
    elseif ind(i)>302 && ind(i)<318
        n(i)=4;
    elseif ind(i)>320 && ind(i)<335
        n(i)=5;
    elseif ind(i)>340 && ind(i)<358
        n(i)=6;
    elseif ind(i)>360 && ind(i)<380
        n(i)=7;
    elseif ind(i)>385 && ind(i)<402
        n(i)=8;
    elseif ind(i)>410 && ind(i)<423
        n(i)=9;
    elseif ind(i)>430 && ind(i)<448
        n(i)=10;
    elseif ind(i)>455 && ind(i)<475
        n(i)=11;
    elseif ind(i)>490 && ind(i)<502
        n(i)=12;
    else
        n(i)=0;
    end
    i=i+1;
end
f6=figure;
scatter(t1,n,'filled');
grid on;
xlabel('Time');
ylabel('notes');
title('Transcription of audio');

%%    end    %%
