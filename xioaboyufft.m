%时频分析
% 原始信号
fs=1000;
f1=50;
f2=100;
t=0:1/fs:1;
s = sin(2*pi*f1*t)+sin(2*pi*f2*t);%+randn(1, length(t));
% s = chirp(t,20,0.3,300);
s = chirp(t,20,1,500,'q');
figure
plot(t, s)
% 连续小波变换时频图
wavename='cmor3-3';
totalscal=256;
Fc=centfrq(wavename); % 小波的中心频率
c=2*Fc*totalscal;
scals=c./(1:totalscal);
f=scal2frq(scals,wavename,1/fs); % 将尺度转换为频率
coefs=cwt(s,scals,wavename); % 求连续小波系数
figure
subplot(221);
imagesc(t,f,abs(coefs));
set(gca,'YDir','normal')
%colorbar;
xlabel('时间 t/s');
ylabel('频率 f/Hz');
title('小波时频图');
% 短时傅里叶变换时频图
subplot(222);
s(s,256,250,256,fs);
% 时频分析工具箱里的短时傅里叶变换
f = 0:fs/2;
stft= stft(s');
stft = stft(1:floor(length(s)/2), :);
figure
subplot(223);
imagesc(t, f, abs(stft));
set(gca,'YDir','normal')
%colorbar;
xlabel('时间 t/s');
ylabel('频率 f/Hz');
title('短时傅里叶变换时频图');
