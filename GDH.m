% 读取源图像
srcimg = imread('airplane.tiff');
srcimg = rgb2gray(srcimg); % 源图像灰度化
[srcHist,srcX] = imhist(srcimg,256);% 计算源图像的直方图

dstimg = imread('airplane.tiff'); % 读取直方图规定化的目标图像
dstimg = rgb2gray(dstimg);% 目标图像灰度化
[dstHist,dstX] = imhist(dstimg,256); % 计算目标图像的灰度直方图

g2 = histeq(srcimg,dstHist);
[g2Hist,g2X] = imhist(g2);
