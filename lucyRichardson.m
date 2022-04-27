%%
clc,clear,close all
% 原始图像
f = checkerboard(8);
% 噪声滤波器
PSF = fspecial('motion', 7, 45);
% 退化图像
gb = imfilter( f, PSF, 'circular' );
% 高斯滤波
noise = imnoise( zeros(size(f)), 'gaussian', 0, 0.001 );
% 将噪声加到原图上
g = gb + noise;

figure
subplot(221), imshow( f ), title('原图')
subplot(222), imshow( noise, [] ), title('高斯噪声')
subplot(223), imshow( gb ), title('退化图像')
subplot(224), imshow( g ), title('退化图像加高斯噪声')

%% 图像复原
% 噪信比默认为0，即信噪比相当于无穷大，即直接逆滤波的过程
fr1 = deconvwnr( g, PSF );

Sn = abs( fft2(noise) ) .^2;
nA = sum( Sn(:) ) / numel( noise );
Sf = abs( fft2(f) ) .^2;
fA = sum( Sf(:) ) / numel( f );
R = nA / fA;  % 平均噪信比计算
% 使用常数比例的维纳滤波进行复原
fr2 = deconvwnr( g, PSF, R );
% 使用自相关函数的维纳滤波进行复原
Ncorr = fftshift( real(ifft2(Sn)) );
Fcorr = fftshift( real(ifft2(Sf)) );
fr3 = deconvwnr( g, PSF, Ncorr, Fcorr );

figure
subplot(221), imshow( g ), title('模糊的运动图像')
subplot(222), imshow( fr1, [] ), title('直接逆滤波复原图像')
subplot(223), imshow( fr2, [] ), title('常数比例维纳滤波复原图像')
subplot(224), imshow( fr3, [] ), title('使用自相关函数的维纳滤波复原图像')

%% 约束最小二乘滤波
% 4 约等于  64*64*0.001
fr1 = deconvreg( g, PSF, 4 );

fr2 = deconvreg( g, PSF, 0.4, [1e-7, 1e7] );

figure
subplot(221), imshow( g ), title('模糊的运动图像')
subplot(222), imshow( fr1, [] ), title('仅噪声功率参数的正则滤波')
subplot(223), imshow( fr2, [] ), title('包含噪声和gamma范围的正则滤波')

%% Lucy-Richardson算法的迭代非线性复原
g = checkerboard(8);
%g=imread('peppers.tiff');
ori = g;
PSF = fspecial( 'gaussian', 7, 10 );
SD = 0.01;
g = imnoise( imfilter(g, PSF), 'gaussian', 0, SD^2 );

damper = 10*SD;
lim = ceil( size(PSF, 1) / 2 );
weight = zeros( size(g) );
weight( lim+1:end-lim, lim+1:end-lim ) = 1;
numit = 5;
f5 = deconvlucy( g, PSF, numit, damper, weight );
numit = 20;
f20 = deconvlucy( g, PSF, numit, damper, weight );
numit = 50;
f50 = deconvlucy( g, PSF, numit, damper, weight );
numit = 100;
f100 = deconvlucy( g, PSF, numit, damper, weight );


figure
subplot(231), imshow( ori ), title('原始图像')
subplot(232), imshow( g, [] ), title('加两次高斯噪声的图')
subplot(233), imshow( f5, [] ), title('LR 5次迭代')
subplot(234), imshow( f20 ), title('LR 20次迭代')
subplot(235), imshow( f50, [] ), title('LR 50次迭代')
subplot(236), imshow( f100,[] ), title('LR 100次迭代')
