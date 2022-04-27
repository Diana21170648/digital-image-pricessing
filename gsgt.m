
imgrgb = imread('airplane.tiff'); %读取彩色图像
f = rgb2gray(imgrgb); %将rgb图像转换成灰度图像
subplot(3, 2, 1);
imshow(f)
title('原始图像')
 
%高斯高通滤波
I = double(f);
g = fft2(I);%二维傅立叶变换
g = fftshift(g);%频移
[M, N] = size(g);
D0 = [10,30,80,150];%截止频率为5
m = fix(M / 2); n = fix(N / 2);
 
for i = 1:M
    for j = 1:N
        D = sqrt((i - m)^2 + (j - n)^2);
        H1 = exp(-(D.^2) ./ (2 * (D0(1)^2)));
        result1(i, j) = (1 - H1) * g(i, j);
    end
end
result1 = ifftshift(result1);
J1 = ifft2(result1);
J2 = uint8(real(J1));
subplot(3, 2, 3);
imshow(J2)
title('高斯高通滤波 D0=10')

for i = 1:M
 
    for j = 1:N
        D = sqrt((i - m)^2 + (j - n)^2);
        H2 = (exp(-(D.^2) ./ (2 * (D0(2)^2))));
        result2(i, j) = (1 - H2) * g(i, j);
    end
end
 
result2 = ifftshift(result2);
J1 = ifft2(result2);
J2 = uint8(real(J1));
subplot(3, 2, 4);
imshow(J2)
title('高斯高通滤波 D0=30')

for i = 1:M
 
    for j = 1:N
        D = sqrt((i - m)^2 + (j - n)^2);
        H3 = exp(-(D.^2) ./ (2 * (D0(3)^2)));
        result3(i, j) = (1 - H3) * g(i, j);
    end
end 
result3 = ifftshift(result3);
J1 = ifft2(result3);
J2 = uint8(real(J1));
subplot(3, 2, 5);
imshow(J2)
title('高斯高通滤波 D0=80')
for i = 1:M
 
    for j = 1:N
        D = sqrt((i - m)^2 + (j - n)^2);
        H4 = exp(-(D.^2) ./ (2 * (D0(4)^2)));
        result4(i, j) = (1 - H4) * g(i, j);
    end
end
 
result4 = ifftshift(result4);
J1 = ifft2(result4);
J2 = uint8(real(J1));
subplot(3, 2, 6);
imshow(J2)
title('高斯高通滤波 D0=150')


