close all;clear all;clc;
I=imread('airplane.tiff');
I=im2double(I);
[M,N]=size(I);
B=zeros(size(I));
for x=2:M-1
    for y=2:N-1
        B(x,y)=I(x+1,y)+I(x-1,y)+I(x,y+1)+I(x,y-1)-4*I(x,y);
    end
end
I=im2uint8(I);
B=im2uint8(B);
subplot(121);imshow(I);title('原图');
subplot(122);imshow(B);title('拉普拉斯算子后的图');
