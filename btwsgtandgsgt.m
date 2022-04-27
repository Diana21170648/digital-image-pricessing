% 1、对lena图像采用高频强调滤波增强方法，
%    并分析方法的效果。（理想、巴特沃斯、高斯）
%    其结果好不好？能否有改善的方法？
%
% @author: jackma
% @time:   2020-10-19 10:15
% @URL:    www.jackrma.com
% @Copyright:博客所有权归本人和CSDN所有，如有转载请在显著位置给出博文
%            链接和作者姓名，否则本人将付诸法律。
% @edit:   
 
% 1. Lena图像采用高频强调滤波增强方法
clc
clear
imgrgb = imread('airplane.tiff'); %读取彩色图像
f = rgb2gray(imgrgb); %将rgb图像转换成灰度图像
subplot(2, 2, 1);
imshow(f)
title('原始图像')
 
%高斯高通滤波
I = double(f);
g = fft2(I);%二维傅立叶变换
g = fftshift(g);%频移
[M, N] = size(g);
D0 = 5;%截止频率为5
m = fix(M / 2); n = fix(N / 2);
 
for i = 1:M
 
    for j = 1:N
        D = sqrt((i - m)^2 + (j - n)^2);
        H = exp(-(D.^2) ./ (2 * (D0^2)));
        result(i, j) = (1 - H) * g(i, j);
    end
 
end
 
result = ifftshift(result);
J1 = ifft2(result);
J2 = uint8(real(J1));
subplot(2, 2, 2);
imshow(J2)
title('高斯高通滤波后的图像')
 
%高频强调滤波
F = 0.5 + 0.75 * (1 - H);
G = F * g;
result2 = ifftshift(G);
J3 = ifft2(result2);
J4 = uint8(real(J3));
subplot(2, 2, 3)
imshow(J4)
title('高频强调滤波后的图像')
 
%对高频强调滤波后图像进行直方图均衡化
J5 = histeq(J4, 256);
J5 = uint8(J5);
subplot(2, 2, 4);
imshow(J5)
title('直方图均衡化后的图像')
 
%%%%%巴特沃斯高通滤波
figure('Name', '图像加入巴特沃斯高通滤波'); %标题
[M, N] = size(f);
a = fft2(f);
a = fftshift(a);
m1 = fix(M / 2); n1 = fix(N / 2);
 
for u = 1:M
 
    for v = 1:N
        D1 = sqrt((u - m1)^2 + (v - n1)^2);
 
        if D1 == 0
            H1(u, v) = 0;
        else
            %    H(u,v)=1/(1+0.414*(500/D1)^4);%截至频率为500
            H1(u, v) = 1 / (1 + (500 / D1)^4); %2阶巴特沃斯高通滤波器，截至频率为500
        end
 
    end
 
end
 
F1 = H1 .* a;
F1 = ifftshift(F1);
I2 = abs(ifft2(F1));
subplot(2, 2, 1);
imshow(f)
title('原始图像')
subplot(2, 2, 2);
imshow(I2)
title('巴特沃斯高通滤波后的图像')
 
%高频强调滤波
FF = 0.5 + 0.75 * (1 - H1);
G1 = FF .* a;
result3 = ifftshift(G1);
J8 = ifft2(result3);
J9 = uint8(real(J8));
subplot(2, 2, 3)
imshow(J9)
title('高频强调滤波后的图像')
 
%对高频强调滤波后图像进行直方图均衡化
J10 = histeq(J9, 256);
J10 = uint8(J10);
subplot(2, 2, 4);
imshow(J10)
title('直方图均衡化后的图像')