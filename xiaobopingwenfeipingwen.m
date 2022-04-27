%%小波平稳与非平稳信号时频图（平稳/非平稳信号）和短时傅里叶变化时频图（平稳/非平稳信号）
%时频分析
% 原始信号
fs=1000;
f1=50;
f2=100;
t=0:1/fs:1;
s = sin(2*pi*f1*t)+sin(2*pi*f2*t);%+randn(1, length(t));
% s = chirp(t,20,0.3,300);
%s = chirp(t,20,1,500,'q');
figure(1)
subplot(211);
plot(t, s)
% 连续小波变换时频图
wavename='cmor3-3';
totalscal=256;
Fc=centfrq(wavename); % 小波的中心频率
c=2*Fc*totalscal;
scals=c./(1:totalscal);
f=scal2frq(scals,wavename,1/fs); % 将尺度转换为频率
coefs=cwt(s,scals,wavename); % 求连续小波系数
figure(2);
subplot(211);
imagesc(t,f,abs(coefs));
set(gca,'YDir','normal')
%colorbar;
xlabel('时间 t/s');
ylabel('频率 f/Hz');
title('平稳小波时频图');
%% 短时傅里叶变换
N=size(t,2);
wlen=500;%设置窗口长度。窗口越长时间分辨率越差，频率分辨率越好。
hop=1;%每次平移的步长，最小为1。越小图像时间精度越好，但计算量大。
x=wkeep1(x,N+1*wlen);%中间截断
h=hamming(wlen);%设置海明窗的窗长
[B, F, T, P] = spectrogram(x,h,wlen-hop,N,fs);   % B是F行T列的频率峰值，P是对应的能量谱密度
figure(4);
subplot(211);
imagesc(T,F,P);
set(gca,'YDir','normal')
colorbar;
xlabel('时间 t/s');
ylabel('频率 f/Hz');
title('平稳短时傅里叶时频图');
toc;

%时频分析
% 原始信号
fs=1000;
f1=50;
f2=100;
t=0:1/fs:1;
s = sin(2*pi*f1*t)+sin(2*pi*f2*t);%+randn(1, length(t));
% s = chirp(t,20,0.3,300);
s = chirp(t,20,1,500,'q');
figure(1);
subplot(212);
plot(t, s)
% 连续小波变换时频图
wavename='cmor3-3';
totalscal=256;
Fc=centfrq(wavename); % 小波的中心频率
c=2*Fc*totalscal;
scals=c./(1:totalscal);
f=scal2frq(scals,wavename,1/fs); % 将尺度转换为频率
coefs=cwt(s,scals,wavename); % 求连续小波系数
figure(2);
subplot(212);
imagesc(t,f,abs(coefs));
set(gca,'YDir','normal')
%colorbar;
xlabel('时间 t/s');
ylabel('频率 f/Hz');
title('非平稳小波时频图');

%% 短时傅里叶变换
N=size(t,2);
wlen=500;%设置窗口长度。窗口越长时间分辨率越差，频率分辨率越好。
hop=1;%每次平移的步长，最小为1。越小图像时间精度越好，但计算量大。
x=wkeep1(x,N+1*wlen);%中间截断
h=hamming(wlen);%设置海明窗的窗长
[B, F, T, P] = spectrogram(x,h,wlen-hop,N,fs);   % B是F行T列的频率峰值，P是对应的能量谱密度
figure(4);
subplot(212);
imagesc(T,F,P);
set(gca,'YDir','normal')
colorbar;
xlabel('时间 t/s');
ylabel('频率 f/Hz');
title('非平稳短时傅里叶时频图');
toc;