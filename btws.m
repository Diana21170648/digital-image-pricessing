%巴特沃斯高通滤波
imgrgb = imread('airplane.tiff'); %读取彩色图像
f = rgb2gray(imgrgb); %将rgb图像转换成灰度图像
subplot(3, 2, 1);
imshow(f)
title('原始图像')

%figure('Name', '图像加入巴特沃斯高通滤波'); %标题
[M, N] = size(f);
a = fft2(f);
a = fftshift(a);
D0 = [10,30,80,150];%截止频率为5
m1 = fix(M / 2); 
n1 = fix(N / 2); 
result1=zeros(M, N);  %预先分配内存空间，提高运行速率
result2=zeros(M, N);  %预先分配内存空间，提高运行速率
result3=zeros(M, N);  %预先分配内存空间，提高运行速率
result4=zeros(M, N);  %预先分配内存空间，提高运行速率

for u = 1:M
    for v = 1:N
        D1 = sqrt((u - m1)^2 + (v - n1)^2);
            H1 = 1 / (1 + ((D0(1) / D1)^4));%2阶巴特沃斯高通滤波器，截至频率为500
             result1(u,v)= g(u,v)*H1; 
    end
end
F1 = H1 .* a;
F1 = ifftshift(F1);
I2 = abs(ifft2(F1));
subplot(3, 2, 3);
imshow(I2)
title('10巴特沃斯高通滤波后的图像')