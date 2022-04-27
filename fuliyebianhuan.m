%平稳信号傅里叶变换  
fs=1000;%采样频率  是最高频率的两倍以上
f1=50;
f2=100;
t=0:1/fs:1;
%t=1/fs;%采样间隔时长
T=2;%采样窗口长度
%A=10;%幅度
N=T*fs% 采样点数 T/t
%x=[0:(1/fs):T];%周期
%len=length(x)-1;%x=N
y=sin(2*pi*f1*t)+sin(2*pi*f2*t);
%y1=sin(2*pi*f1*t)+sin(2*pi*f2*t);
%y = chirp(t,20,1,500,'q');
figure
subplot(211);
%subplot(3,1,1);
%plot(t,y);
%xlabel('时间/s');
%ylabel('幅度');
%subplot(3,1,2);
[f,y]=get_fft(y,fs,N);
plot(f,y);
xlabel('f/Hz');
ylabel('幅度');

%非平稳信号傅里叶变换 
fs=1000;%采样频率  是最高频率的两倍以上
f1=50;
f2=100;
t=0:1/fs:1;
%t=1/fs;%采样间隔时长
T=2;%采样窗口长度
%A=10;%幅度
N=T*fs% 采样点数 T/t
%x=[0:(1/fs):T];%周期
%len=length(x)-1;%x=N
y=sin(2*pi*f1*t)+sin(2*pi*f2*t);
%y1=sin(2*pi*f1*t)+sin(2*pi*f2*t);
y = chirp(t,20,1,500,'q');

subplot(212);
%subplot(3,1,1);
%plot(t,y);
%xlabel('时间/s');
%ylabel('幅度');
%subplot(3,1,2);
[f,y]=get_fft(y,fs,N);
plot(f,y);
xlabel('f/Hz');
ylabel('幅度');
function [f, spectrum ] = get_fft(s,Fs,L)
%GAN_FFT 此处显示有关此函数的摘要
%   此处显示详细说明
y=fft(s);
p2=abs(y/L);
p1=p2(1:L/2+1);
p1(2:end-1)=2*p1(2:end-1);
f = Fs*(0:(L/2))/L;
spectrum=p1;
 
end
%plot(f,spectrum);


