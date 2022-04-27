%理想低通
I = imread('airplane.tiff');
I=rgb2gray(I);
figure(1);
subplot(321),imshow(I);
title('原图像');
%I=imnoise(I,'gaussian');%%加入高斯白噪声
%subplot(222),imshow(I);
%title('加入噪声后的图像');
s=fftshift(fft2(I));
%subplot(223), imshow(log(abs(s)),[]); 
%title('图像傅里叶变换取对数所得频谱');
[a,b]=size(s);
a0=round(a/2);
b0=round(b/2);
d=[10,30,80,150];
for i=1:a 
    for j=1:b 
        distance=sqrt((i-a0)^2+(j-b0)^2);
        if distance<=d(1)
            h=1;
        else
            h=0;
        end
        s1(i,j)=h*s(i,j);
    end
end
output1=uint8(real(ifft2(ifftshift(s1))));
subplot(323),imshow(output1);
title('10理想低通滤波所得图像'); 

for i=1:a 
    for j=1:b 
        distance=sqrt((i-a0)^2+(j-b0)^2);
        if distance<=d(2)
            h=1;
        else
            h=0;
        end
        s2(i,j)=h*s(i,j);
    end
end
output2=uint8(real(ifft2(ifftshift(s2))));
subplot(324),imshow(output2);
title('30理想低通滤波所得图像'); 

for i=1:a 
    for j=1:b 
        distance=sqrt((i-a0)^2+(j-b0)^2);
        if distance<=d(3)
            h=1;
        else
            h=0;
        end
        s3(i,j)=h*s(i,j);
    end
end
output3=uint8(real(ifft2(ifftshift(s3))));
subplot(325),imshow(output3);
title('80理想低通滤波所得图像'); 

for i=1:a 
    for j=1:b 
        distance=sqrt((i-a0)^2+(j-b0)^2);
        if distance<=d(4)
            h=1;
        else
            h=0;
        end
        s4(i,j)=h*s(i,j);
    end
end
output4=uint8(real(ifft2(ifftshift(s4))));
subplot(326),imshow(output4);
title('150理想低通滤波所得图像'); 

