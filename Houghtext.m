%用Hough变换提取直线
%在对图形进行直线检测之前，需要对图形进行边缘检测、二值化处理，用拉普拉斯算法或canny提取边缘
%Hough变换思想
%1．在ρ、θ的极值范围内对其分别进行 m，n 等分，设一个二维数组的下标与 ρi 、 θj 的取值 对应；
%2．对图像上的所有边缘点作 Hough 变换，求每个点在 θj (j＝0，1，…，n)Hough 变换后的 ρi，判断( ρi 、 θj )与哪个数组元素对应，则让该数组元素值加 1； 
%3．比较数组元素值的大小，最大值所对应的( ρi 、 θj )就是这些共线点对应的直线方程的参数。 
close all;clear;clc;
a=imread('车道检测.jpg');
a=rgb2gray(a);
figure;
a=a(:,:,1);
subplot(221);%原图
imshow(a); 

bw1=LapuLas(a);%调用自己编写的拉普拉斯算法进行边缘检测
subplot(222);%拉普拉斯检测结果
imshow(bw1);
%bw1=edge(a,'canny',0.2);
%imshow(bw1);

%subplot (223);%Hough空间
[H,theta,rho]=naiveHough(bw1);%调用自己编写的Hough函数,H为Hough变换矩阵，theta和rho为霍夫变换的角度和半径
imshow(H,[],'XData',theta,'YData',rho, 'InitialMagnification','fit');%画出Hough空间
xlabel('\theta (degrees)'), ylabel('\rho');
axis on, axis square, hold on;
P = houghpeaks(H,6);%提取6个极值点
x = theta(P(:,2));
y = rho(P(:,1));
plot(x,y,'s','color','red');%标出极值点
title('Hough空间');%hough变换的结果

lines = houghlines(bw1,theta,rho,P);%提取线段
subplot(224)
imshow(a),hold on;%为了在原图上画出直线段
for k = 1:length(lines)
xy = [lines(k).point1; lines(k).point2];
 plot(xy(:,1),xy(:,2),'LineWidth',1,'Color','green');%画出线段
plot(xy(1,1),xy(1,2),'o','LineWidth',1,'Color','yellow');%起点
plot(xy(2,1),xy(2,2),'o','LineWidth',1,'Color','red');%终点
end
%Hough变换函数
function [ Hough, theta_range, rho_range ] = naiveHough(I)
[rows, cols] = size(I);
theta_maximum = 90;
rho_maximum = floor(sqrt(rows^2 + cols^2)) - 1;
theta_range = -theta_maximum:theta_maximum - 1;
rho_range = -rho_maximum:rho_maximum;
Hough = zeros(length(rho_range), length(theta_range));
for row = 1:rows
    for col = 1:cols
        if I(row, col) > 0 %only find: pixel > 0
            x = col - 1;
            y = row - 1;
            for theta = theta_range
                rho = round((x * cosd(theta)) + (y * sind(theta)));  
                rho_index = rho + rho_maximum + 1;
                theta_index = theta + theta_maximum + 1;
                Hough(rho_index, theta_index) = Hough(rho_index, theta_index) + 1;
            end
        end
    end
end
end
%拉普拉斯边缘提取
function [p] = LapuLas(e)
r=e;
[m,n]=size(e);
%图像二值化后，拉普拉斯算法提取边缘
dd=sum(sum(e)/(m*n));
for x=1:m
    for y=1:n
        if(r(x,y)>=dd)
            r(x,y)=255;
        else
            r(x,y)=0;
        end
    end
end
subplot(2,2,3),imshow(r);title('二值化后图像');
p=zeros(m-1,n-1);
for ii=2:m-1
    for jj=2:n-1
        p(ii,jj)=r(ii,jj+1)+r(ii,jj-1)+r(ii+1,jj)+r(ii-1,jj)-4*(r(ii,jj));
    end
end
p=uint8(p);
end



%R-Hough变换 BW5 --------------------------------------------------------------------
%function Untitled_13_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%P-Hough变换 BW6 --------------------------------------------------------------------
%function Untitled_14_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% S-Hough变换 BW7--------------------------------------------------------------------
%function Untitled_15_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%L-Hough变换 BW1--------------------------------------------------------------------
%function Untitled_16_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
