clc;
clear;
close all;

% 读取源图像
p = imread('airplane.tiff');
%srcimg = rgb2gray(srcimg); % 源图像灰度化
[srcHist,srcX] = imhist(p,256);% 计算源图像的直方图
srcHist = srcHist/numel(p);% 直方图归一化
cps = zeros(256,1,'double');% 计算源图像的灰度累积概率
for i=1:1:256
    cps(i) = sum(srcHist(1:i));
end
cps = 255*cps; % 构建源图像灰度均衡化的变换函数
cps = uint8(cps);

P1 = imread('peppers.tiff'); % 读取直方图规定化的目标图像
%dstimg = rgb2gray(dstimg);% 目标图像灰度化
[dstHist,dstX] = imhist(P1,256); % 计算目标图像的灰度直方图
dstHist = dstHist/sum(dstHist); % 灰度直方图归一化
cpd = zeros(256,1,'double');% 计算目标图像的灰度分布累积概率
for i=1:1:256
    cpd(i) = sum(dstHist(1:i));
end
cpd = cpd*255; % 构建目标图像的灰度均衡化变换函数
cpd = uint8(cpd);

srcl=zeros(256,1,'uint8');
minv = 256;
for i = 1:1:256 
    minv =256;
    for j = 1:1:256
        if minv > abs(cps(i)-cpd(j))
            minv =  abs(cps(i)-cpd(j));
            srcl(i) =j;
        end
    end
end
[width,height] = size( p);
gray1 = p;
for i=1:1:width
    for j = 1:1:height
        gray1(i,j)=srcl(p(i,j)+1);
    end
end

[g1Hist,g1X] = imhist(gray1);
% g1Hist = g1Hist/sum(g1Hist);
g2 = histeq(p,dstHist);
[g2Hist,g2X] = imhist(g2);
% g2Hist = g2Hist/sum(g2Hist);
figure(1),
subplot(2,4,1),imshow(p),title('原图');
subplot(2,4,2),imshow(P1),title('匹配图');
subplot(2,4,3),imshow(gray1),title('结果1');
subplot(2,4,4),imshow(g2),title('结果2');
subplot(2,4,5),stem(srcX,srcHist),title('原图直方图');
subplot(2,4,6),stem(dstX,dstHist),title('匹配图直方图');
subplot(2,4,7),stem(g1X,g1Hist),title('结果1直方图');
subplot(2,4,8),stem(g2X,g2Hist),title('结果2直方图');
