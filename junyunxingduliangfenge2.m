%均匀性度量分割
clc;
clear all;
close all;
I=imread('cameraman.tif');
subplot(121),imshow(I);title('Original');
% subplot(222),imhist(Image);title(' histogram');
%Initial threshold
I=double(I);
minValue=min(min(I));
maxValue=max(max(I));
[row,col]=size(I);
Th=minValue+1; %给定初始阈值
perfactValue=10000000000; %假设初始为无穷大

for m=minValue+1:maxValue-1
k1=1;k2=1;
for i=1:row
for j=1:col

if I(i,j)<m
C1(1,k1)=I(i,j);k1=k1+1; %C1类
else
C2(1,k2)=I(i,j);k2=k2+1; %C2类
end

end
end

%对C1类求均值，方差，分布概率
average1=mean(C1); %均值1
variance1=0;
for i=1:k1-1
variance1=variance1+(C1(1,i)-average1)^2; %C1类的方差
end
variance1 = variance1/(k1-1);

p1=(k1-1)/(row*col); %C1类的分布概率

%对C2类求均值，方差，分布概率
average2=mean(C2); %均值2
variance2=0;
for i=1:k2-1
variance2=variance2+(C2(1,i)-average2)^2; %C2类的方差
end
p2=(k2-1)/(row*col); %C2类的分布概率
variance2 = variance2/(k2-1);

newValue=p1*variance1+p2*variance2;
if (newValue<perfactValue)
Th=m;
perfactValue=newValue;
end
end 
% Th=82;
for i=1:row
for j=1:col
if I(i,j) >= Th
G(i,j)=255;
else
G(i,j)=0;
end
end
end
subplot(122),imshow(G);title('segmentation');