%有约束最小二乘滤波复原
I = im2double(imread('peppers.tiff'));
[hei,wid,~] = size(I);
subplot(2,2,1),imshow(I);
title('原图像');
% 模拟运动模糊.
LEN = 21;
THETA = 11;
PSF = fspecial('motion', LEN, THETA);%产生运动模糊算子，即点扩展函数
blurred = imfilter(I, PSF, 'conv', 'circular');
subplot(2,2,2), imshow(blurred); title('模糊图像');
Pf = psf2otf(PSF,[hei,wid]);%退化函数的FFT
% 添加加性噪声
noise_mean = 0;
noise_var = 0.00001;
blurred_noisy = imnoise(blurred, 'gaussian',noise_mean, noise_var);
subplot(2,2,3), imshow(blurred_noisy)
title('带运动模糊和噪声图像')
p = [0 -1 0;-1 4 -1;0 -1 0];%拉普拉斯模板
P = psf2otf(p,[hei,wid]);
gama = 0.001;
If = fft2(blurred_noisy);
numerator = conj(Pf);%conj函数，用于求一个复数的复共轭
denominator = Pf.^2 + gama*(P.^2);
deblurred2 = ifft2( numerator.*If./ denominator );%约束最小二乘方滤波在频率域中的表达式
subplot(2,2,4), imshow(deblurred2)
title('约束最小二乘方滤波后图像');
