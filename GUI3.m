
function varargout = GUI3(varargin)
% GUI3 MATLAB code for GUI3.fig
%      GUI3, by itself, creates a new GUI3 or raises the existing
%      singleton*.
%
%      H = GUI3 returns the handle to a new GUI3 or the handle to
%      the existing singleton*.
%
%      GUI3('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI3.M with the given input arguments.
%
%      GUI3('Property','Value',...) creates a new GUI3 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI3_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI3_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI3

% Last Modified by GUIDE v2.5 19-Apr-2022 20:51:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI3_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI3_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GUI3 is made visible.
function GUI3_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI3 (see VARARGIN)

% Choose default command line output for GUI3
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI3 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = GUI3_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%读入图像并在坐标里面显示 --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname,filterindex]=...
uigetfile({'*.*';'*.bmp';'*.tif';'*.png';'*.jpg';'*.jpeg'},'select picture');
str=[pathname filename];  
s=str;
handles.filebig=filterindex;
if filterindex==0
%msgbox('选择图像失败！','error');
return
else   
im1=imread(str);  
end 
axes(handles.axes1);  
imshow(im1);  %显示图片
handles.img=im1;
guidata(hObject,handles);


% 图像分割目录--------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%边缘检测目录 --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Hough变换(这利用拉普拉斯算子）--------------------------------------------------------------------
function Untitled_4_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%梯度算子目录 --------------------------------------------------------------------
function Untitled_5_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% 拉普拉斯算子--------------------------------------------------------------------
function Untitled_6_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global BW4
I = getimage(gca);
%I=rgb2gray(I);
%I=imread('airplane.tiff');
I=im2double(I);
[M,N]=size(I);
BW1=zeros(size(I));
for x=2:M-1
    for y=2:N-1
        BW1(x,y)=I(x+1,y)+I(x-1,y)+I(x,y+1)+I(x,y-1)-4*I(x,y);
    end
end
%I=im2uint8(I);
BW1=im2uint8(BW1);
%BW = edge(I,'laplacian'); %利用canny算子进行边缘检测
axes(handles.axes2);  
imshow(BW1);  %显示图片
handles.img=BW1;
guidata(hObject,handles);

% LOG算子--------------------------------------------------------------------
function Untitled_7_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = getimage(gca);
I=rgb2gray(I);
BW2 = edge(I,'log'); %利用canny算子进行边缘检测
axes(handles.axes2);  
imshow(BW2);  %显示图片
handles.img=BW2;
guidata(hObject,handles);


% DOG算子(两个log相减）--------------------------------------------------------------------
function Untitled_8_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%I = getimage(gca);
%I=double(I)/255
%I=rgb2gray(I);
%<span style="font-size:18px;">;
I = double(getimage(gca))/255;
gray=rgb2gray(I);
 
sigma1=0.1;
sigma2=0.8;
window=7;
H1=fspecial('gaussian', window, sigma1);
H2=fspecial('gaussian', window, sigma2);
 
DiffGauss=H1-H2;
BW3=imfilter(gray,DiffGauss,'replicate');
BW3=mat2gray(BW3);
 
%figure;imshow(BW3);
%</span>
%BW3 = edge(I,'dog'); %利用canny算子进行边缘检测
axes(handles.axes2);  
imshow(BW3);  %显示图片
handles.img=BW3;
guidata(hObject,handles);

% Canny算子--------------------------------------------------------------------
function Untitled_9_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = getimage(gca);
I=rgb2gray(I);
BW4 = edge(I,'canny'); %利用canny算子进行边缘检测
axes(handles.axes2);  
imshow(BW4);  %显示图片
handles.img=BW4;
guidata(hObject,handles);


% 清除检测结果和所选图像，只清除检测结果会导致数组长度不符合
%--- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1); %指定需要清空的坐标轴
cla reset;
axes(handles.axes2); %指定需要清空的坐标轴
cla reset;



% roberts梯度算子--------------------------------------------------------------------
function Untitled_10_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = getimage(gca);
I=rgb2gray(I);
BW5 = edge(I,'roberts'); %利用canny算子进行边缘检测
axes(handles.axes2);  
imshow(BW5);  %显示图片
handles.img=BW5;
guidata(hObject,handles);

%prewitt梯度算子 --------------------------------------------------------------------
function Untitled_11_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = getimage(gca);
I=rgb2gray(I);
BW6 = edge(I,'prewitt'); %利用canny算子进行边缘检测
axes(handles.axes2);  
imshow(BW6);  %显示图片
handles.img=BW6;
guidata(hObject,handles);

% sobel梯度算子--------------------------------------------------------------------
function Untitled_12_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = getimage(gca);
I=rgb2gray(I);
BW7 = edge(I,'sobel'); %利用canny算子进行边缘检测
axes(handles.axes2);  
imshow(BW7);  %显示图片
handles.img=BW7;
guidata(hObject,handles);


% P-Hough--------------------------------------------------------------------
function Untitled_13_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = getimage(gca);
I=rgb2gray(I);             %rgb转灰度图
%subplot(2,2,1);               %显示灰度图
imshow(I);
%title('灰度图像');
axis on;
BW=edge(I,'prewitt');
%subplot(2,2,2);
%imshow(BW);
%title('prewitt算子边缘检测后图像');
axis on;
[H,T,R]=hough(BW);
%subplot(2,2,3);
%imshow(H,[],'XData',T,'YData',R,'InitialMagnification','fit');
%title('霍夫变换图');
xlabel('\theta'),ylabel('\rho');
axis on , axis normal, hold on;
P=houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
x=T(P(:,2));y=R(P(:,1));
plot(x,y,'s','color','white');
lines=houghlines(BW,T,R,P,'FillGap',5,'MinLength',7);
axes(handles.axes2);  
imshow(I);  %显示图片
handles.img=I;
guidata(hObject,handles);
%subplot(2,2,4);imshow(rotI);
title('霍夫变换图像检测');
axis on;
hold on;
max_len=0;
for k=1:length(lines)
xy=[lines(k).point1;lines(k).point2];
plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
len=norm(lines(k).point1-lines(k).point2);
if(len>max_len)
max_len=len;
xy_long=xy;
end
end
plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','cyan');


%L-Hough --------------------------------------------------------------------
function Untitled_14_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a = getimage(gca);
a=rgb2gray(a);
%figure;
%a=a(:,:,1);
%subplot(221);%原图
%imshow(a); 

bw1=LapuLas(a);%调用自己编写的拉普拉斯算法进行边缘检测
%subplot(222);%拉普拉斯检测结果
%imshow(bw1);
%bw1=edge(a,'canny',0.2);
%imshow(bw1);

%subplot (223);%Hough空间
[H,theta,rho]=naiveHough(bw1);%调用自己编写的Hough函数,H为Hough变换矩阵，theta和rho为霍夫变换的角度和半径
%imshow(H,[],'XData',theta,'YData',rho, 'InitialMagnification','fit');%画出Hough空间
%xlabel('\theta (degrees)'), ylabel('\rho');
%axis on, axis square, hold on;
P = houghpeaks(H,6);%提取6个极值点
%x = theta(P(:,2));
%y = rho(P(:,1));
%plot(x,y,'s','color','red');%标出极值点
%title('Hough空间');%hough变换的结果

lines = houghlines(bw1,theta,rho,P);%提取线段
axes(handles.axes2);  
imshow(a);  %显示图片
handles.img=a;
guidata(hObject,handles);
%subplot(224)
%imshow(a),
hold on;%为了在原图上画出直线段
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
%subplot(2,2,3),imshow(r);title('二值化后图像');
p=zeros(m-1,n-1);
for ii=2:m-1
    for jj=2:n-1
        p(ii,jj)=r(ii,jj+1)+r(ii,jj-1)+r(ii+1,jj)+r(ii-1,jj)-4*(r(ii,jj));
    end
end
p=uint8(p);


%图像分割算法菜单 --------------------------------------------------------------------
function Untitled_15_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%区域生长法分割 --------------------------------------------------------------------
function Untitled_16_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = getimage(gca);
seed=[100,220];%选择起始位置
thresh=16;%相似性选择阈值
A=rgb2gray(I);
A=imadjust(A,[min(min(double(A)))/255,max(max(double(A)))/255],[]);
A=double(A);
B=A;
[r,c]=size(B);%图像尺寸，r为行数，c为列数
n=r*c;%计算图像中所包含点的个数
pixel_seed=A(seed(1),seed(2));%原图起始点位置
q=[seed(1) seed(2)];%用q来装载起始位置
top=1;%循环判断flag
M=zeros(r,c);%建立一个与原图形同等大小的矩阵
M(seed(1),seed(2))=1;
count=1;%计数器
while top~=0 %循环结束条件
    r1=q(1,1);
    c1=q(1,2);
    p=A(r1,c1);
    dge=0;
    for i=-1:1
        for j=-1:1
            if r1+i<=r&&r1+i>0&&c1+j<=c&&c1+j>0
                if abs(A(r1+i,c1+j)-p)<=thresh&&M(r1+i,c1+j)~=1
                    top=top+1;
                    q(top,:)=[r1+i c1+j];
                    M(r1+i,c1+j)=1;
                    count=count+1;
                    B(r1+i,c1+j)=1;%满足判定条件将B中对应的点赋为1
                end
                if M(r1+i,c1+j)==0
                    dge=1;%将dge赋为1
                end
            else
                dge=1;%点在图像外，将dge赋为1
            end
        end
    end
    if dge~=1
        B(r1,c1)=A(seed(1),seed(2));%将原图像起始位置灰度值赋予B
    end
    if count >=n
        top =1;
    end
    q=q(2:top,:);
    top=top-1;
end
%subplot(121)
%imshow(A,[])
%title('灰度图像')
%subplot(122)
%imshow(B,[])
%title('区域生长法分割图像')
axes(handles.axes2);  
imshow(B,[]);  %显示图片
handles.img=B,[];
guidata(hObject,handles);

%区域分裂与合并法分割--------------------------------------------------------------------
function Untitled_17_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I1 = getimage(gca);
%I1=rgb2gray(I1);
%I1=imread('rice.png');
S=qtdecomp(I1,0.25);%其中0.25为每个方块所需要达到的最小差值
I2=full(S);
%subplot(121);
%imshow(I1);
%title('原图像')
%subplot(122)
%imshow(I2)
%title('四叉树分解的图像')

axes(handles.axes2);  
imshow(I2);  %显示图片
handles.img=I2;
guidata(hObject,handles);

%迭代阈值法分割--------------------------------------------------------------------
function Untitled_18_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%I=imread('peppers.tiff');
I = getimage(gca);
%I=rgb2gray(I);
I=im2double(I);
level=0.1;
max_I=max(I(:));
min_I=min(I(:));
T1=1/2*(max_I+min_I);
[M,N]=size(I);
A_number=0;
B_number=0;
A_all=0;
B_all=0;
for i=1:M
    for j=1:N
        if(I(i,j)>=T1)
            A_number=A_number+1;
            A_all=A_all+I(i,j);
        else if(I(i,j)<T1)
             B_number=B_number+1;
             B_all=B_all+I(i,j);
            end
        end
    end
    A_ave= A_all/A_number;
    B_ave= B_all/B_number;
    T2=1/2*( A_ave+ B_ave );%*不能省略
end
T2*255
J=im2bw(I,T2);
%figure(1);
%subplot(121);imshow(I);
%subplot(122);imshow(J);
axes(handles.axes2);  
imshow(J);  %显示图片
handles.img=J;
guidata(hObject,handles);

%均匀性度量法分割 （选灰色花朵或者相机的图片）--------------------------------------------------------------------
function Untitled_19_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see G%IDATA)
%I=imread('cameraman.tif');
I = getimage(gca);
%subplot(121),imshow(I);title('Original');
% subplot(222),imhist(Image);title(' histogram');
%Initial threshold
%I=double(I);
I=rgb2gray(I);
minValue=min(min(I));
maxValue=max(max(I));
[row,col]=size(I);
Th=minValue+1; %给定初始阈值
perfactValue=10000000000; %假设初始为无穷大

for m=1:254 %m=minValue+1:maxValue-1
    k1=1;
    k2=1;
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
%subplot(122),imshow(G);title('segmentation');
axes(handles.axes2);  
imshow(G);  %显示图片
handles.img=G;
guidata(hObject,handles);

%聚类法分割 --------------------------------------------------------------------
function Untitled_20_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%im = imread('peppers.tiff');
im = getimage(gca);
    %imshow(im), title('Imput image');  %%转换图像的颜色空间得到样本
    cform = makecform('srgb2lab');
    lab = applycform(im,cform);
    ab = double(lab(:,:,2:3));
    nrows = size(lab,1); ncols = size(lab,2);
    X = reshape(ab,nrows*ncols,2)';
    %figure, scatter(X(1,:)',X(2,:)',3,'filled'),title('image 2');  box on; %显示颜色空间转换后的二维样本空间分布
    %% 对样本空间进行Kmeans聚类
    k = 5; % 聚类个数
    max_iter = 100; %最大迭代次数
    [centroids, labels] = run_kmeans(X, k, max_iter); 

    %% 显示聚类分割结果
    %figure, scatter(X(1,:)',X(2,:)'3,(labels,'filled'),title('image 3'); %显示二维样本空间聚类效果
    hold on; scatter(centroids(1,:),centroids(2,:), 60,'r','filled')
    hold on; scatter(centroids(1,:),centroids(2,:),30,'g','filled')
    box on; hold off;
    %print -dpdf 2D2.pdf

    pixel_labels = reshape(labels,nrows,ncols);
    rgb_labels = label2rgb(pixel_labels);
    %figure, imshow(rgb_labels), title('Segmented Image');
    %print -dpdf Seg.pdf
    axes(handles.axes2);  
imshow(rgb_labels);  %显示图片
handles.img=rgb_labels;
guidata(hObject,handles);
    function [centroids, labels] = run_kmeans(X, k, max_iter)
    % 该函数实现Kmeans聚类
    centroids = X(:,1+round(rand*(size(X,2)-1)));
    labels = ones(1,size(X,2));
    for i = 2:k
        D = X-centroids(:,labels);
        D = cumsum(sqrt(dot(D,D,1)));
        if D(end) == 0, centroids(:,i:k) = X(:,ones(1,k-i+1)); return; end
        centroids(:,i) = X(:,find(rand < D/D(end),1));
        [~,labels] = max(bsxfun(@minus,2*real(centroids'*X),dot(centroids,centroids,1).'));
      end

    %% 标准Kmeans算法
     for iter = 1:max_iter
        for i = 1:k, l = labels==i; centroids(:,i) = sum(X(:,l),2)/sum(l); end
        [~,labels] = max(bsxfun(@minus,2*real(centroids'*X),dot(centroids,centroids,1).'),[],1);
    end
 

    
%大津法分割--------------------------------------------------------------------
function Untitled_21_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%I=imread('peppers.tiff');
I = getimage(gca);
p=zeros(1,256);
for k=1:numel(I)%概率密度
    p(I(k)+1)=p(I(k)+1)+1;
end
p=p/numel(I);
cump=zeros(1,256);%分配运行空间
cump=cumsum(p);
temp=0:255;
temp=temp.*p;
m=cumsum(temp);

m_global=m(256);

var_b=zeros(1,256);
for k=1:256
    if p(k)~=0
        var_b(k)=(m_global*cump(k)-m(k)^2/(cump(k)*(1-cump(k))));%类间方差
    end
end
tempk=0;
countk=0;
for k=1:256
    if var_b(k)==max(var_b)
        tempk=tempk+k;
        countk=countk+1;
    end
end
threshold=round(tempk/countk)-1;
seg_I=zeros(size(I));
for k=1:numel(I)
    if I(k)>=threshold
        seg_I(k)=1;
    end
end
%figure('name','Ostu','NumberTitle','off');
%imshow(seg_I);
axes(handles.axes2);  
imshow(seg_I);  %显示图片
