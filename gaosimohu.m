clc, clear, close all
img_origin=imread('peppers.tiff');%读取图片文件
img_gray=rgb2gray(img_origin);%灰度化
Sigma=300;%σ为高斯模糊半径，半径越大模糊程度越大
for x = 1: 3  % 垂直方向
    for y = 1:3  % 水平方向
        WeightMatrix(x, y)=exp(-((x-1)^2+(y-1)^2)/(2*Sigma^2))/(2*pi*Sigma^2);
    end
end
WeightMatrix=WeightMatrix./sum(sum(WeightMatrix)) %使该3*3矩阵之和等于1
[row, col] = size( img_gray );
for i = 1: row  % 垂直方向
    for j = 1:col  % 水平方向 
        if i==1 || j==1 || i==row || j==col
            img_undist(i, j)=img_gray(i, j);%边缘未处理
        else
            miniMatrix=single(img_gray(i-1:i+1, j-1:j+1));
            img_undist(i, j)=sum(sum( miniMatrix.*WeightMatrix ));%高斯模糊
        end
    end
end
figure(1);
subplot(121);
imshow(img_gray);%显示灰度图
subplot(122);
imshow(img_undist);%显示灰度模糊图