%短时傅里叶变化时频图（平稳/非平稳信号）
tic;
%% 产生信号
fs = 1000;  % 采样频率1KHz
t = 0:1/fs:1-1/fs;
N=size(t,2);
f1=50;
f2=100;
x=sin(2*pi*f1*t)+sin(2*pi*f2*t);
figure(1);
subplot(211);
plot(t,x);
%% 短时傅里叶变换
wlen=500;%设置窗口长度。窗口越长时间分辨率越差，频率分辨率越好。
hop=1;%每次平移的步长，最小为1。越小图像时间精度越好，但计算量大。
x=wkeep1(x,N+1*wlen);%中间截断
h=hamming(wlen);%设置海明窗的窗长
[B, F, T, P] = spectrogram(x,h,wlen-hop,N,fs);   % B是F行T列的频率峰值，P是对应的能量谱密度
figure(2);
subplot(211);
imagesc(T,F,P);
set(gca,'YDir','normal')
colorbar;
xlabel('时间 t/s');
ylabel('频率 f/Hz');
title('平稳短时傅里叶时频图');
toc;


fs = 1000;  % 采样频率1KHz
t = 0:1/fs:1-1/fs;
N=size(t,2);
f1=50;
f2=100;
x=sin(2*pi*f1*t)+sin(2*pi*f2*t);
x = chirp(t,20,1,500,'q');
figure(1);
subplot(212);
plot(t,x);
%% 短时傅里叶变换
wlen=500;%设置窗口长度。窗口越长时间分辨率越差，频率分辨率越好。
hop=1;%每次平移的步长，最小为1。越小图像时间精度越好，但计算量大。
x=wkeep1(x,N+1*wlen);%中间截断
h=hamming(wlen);%设置海明窗的窗长
[B, F, T, P] = spectrogram(x,h,wlen-hop,N,fs);   % B是F行T列的频率峰值，P是对应的能量谱密度
figure(2);
subplot(212);
imagesc(T,F,P);
set(gca,'YDir','normal')
colorbar;
xlabel('时间 t/s');
ylabel('频率 f/Hz');
title('非平稳短时傅里叶时频图');
toc;