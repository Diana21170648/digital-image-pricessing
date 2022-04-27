clc;
clear;
close all;

img=imread('peppers.tiff');
%img=ind2gray(I);
figure,imshow(img),title('原图');
psf=fspecial('disk',8);
res1=deconvblind(img,psf);
figure,imshow(res1),title('盲去卷积10次');
res2=deconvblind(img,psf,20);
figure,imshow(res2),title('盲去卷积20次');
res3=deconvblind(img,psf,30);
figure,imshow(res3),title('盲去卷积30次');
res4=deconvblind(img,psf,50);
figure,imshow(res3),title('盲去卷积50次');
