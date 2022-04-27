%<span style="font-size:18px;">clear;
img = double(imread('车道检测.jpg'))/255;
gray=rgb2gray(img);
 
sigma1=0.1;
sigma2=0.8;
window=7;
H1=fspecial('gaussian', window, sigma1);
H2=fspecial('gaussian', window, sigma2);
 
DiffGauss=H1-H2;
out=imfilter(gray,DiffGauss,'replicate');
out=mat2gray(out);
 
figure;imshow(out);
%</span>