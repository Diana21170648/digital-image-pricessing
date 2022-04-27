I=imread('peppers.tiff');
r=10;%散焦半径r
PSF=fspecial('disk',r);   %得到点扩散函数
I1=imfilter(I,PSF,'symmetric','conv');  %实现散焦模糊

%利用拉普拉斯算子对散焦模糊图像进行二阶微分
I1=double(I1);
h=[1 1 1;1 -8 1;1 1 1];
I2=filter2(h,I1);
%对微分图I2进行自相关计算
R=xcorr2(I2);
R=R/max(R(:));
figure,surfc(R); 