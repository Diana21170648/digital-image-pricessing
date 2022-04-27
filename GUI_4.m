function varargout = GUI_4(varargin)
% GUI_4 MATLAB code for GUI_4.fig
%      GUI_4, by itself, creates a new GUI_4 or raises the existing
%      singleton*.
%
%      H = GUI_4 returns the handle to a new GUI_4 or the handle to
%      the existing singleton*.
%
%      GUI_4('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_4.M with the given input arguments.
%
%      GUI_4('Property','Value',...) creates a new GUI_4 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_4_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_4_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_4

% Last Modified by GUIDE v2.5 19-Apr-2022 21:19:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_4_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_4_OutputFcn, ...
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


% --- Executes just before GUI_4 is made visible.
function GUI_4_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_4 (see VARARGIN)

% Choose default command line output for GUI_4
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_4 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_4_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% 图像复原菜单栏按钮--------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% 逆滤波复原--------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = getimage(gca);
%I=rgb2gray(I); 
x1=I(:,:,1);
x1=double(x1); 
[r,r1]=size(x1); 
y1=fftshift(fft2(x1));
[r,r1]=size(y1); 
%figure(1);
%imshow(x1,[]);
%title('原始的图像') ;
%figure(2);
%imshow(abs(y1),[0,250000]);
%title('原始的图像频谱');
m=1:r; 
m1=1:r1; 
[m,m1]=meshgrid(m,m1);%生成网格空间 
noise=20.*imnoise(zeros(r,r1),'gaussian',0,0.008);%高斯噪声 
%figure(3);  
%subplot(1,2,1);
%imshow(noise,[]);
%title('白噪声') ;
a=double(21/100);%x方向的最大移动量为ra的0.21倍,可调 
b=double(21/100);%y方向的最大移动量为ca的0.21倍,可调 
t=double(88/100);%移动到最大所需的时间默认为0.88
f=ones(r,r1);   
g=(m-r/2-1).*a+(m1-r1/2-1).*b+eps;   
f=t.*sin(pi.*g).*exp(-j.*pi.*g)./(pi.*g);
h=f'.*y1; 
tu=ifft2(h);  
tu=abs(tu)+noise; 
%subplot(1,2,2);
%imshow(tu,[]);
%title('加入白噪声退化的图像')%原图傅立叶变换估计值 
axes(handles.axes2);  
imshow(tu,[]);  %显示图片
handles.img=tu;
guidata(hObject,handles);


y1=h./f'; 
%figure(4)
%subplot(1,2,1);
%imshow(abs(ifft2(y1)),[]);
%title('逆滤波的结果');
axes(handles.axes3);  
imshow(abs(ifft2(y1)),[]);  %显示图片
handles.img=abs(ifft2(y1));
guidata(hObject,handles);



%维纳滤波复原 --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = getimage(gca);
%I=rgb2gray((I); 
x1=I(:,:,1);
x1=double(x1); 
[r,r1]=size(x1); 
y1=fftshift(fft2(x1));
[r,r1]=size(y1); 
%figure(1);
%imshow(x1,[]);
%title('原始的图像') ;
%figure(2);
%imshow(abs(y1),[0,250000]);
%title('原始的图像频谱');
m=1:r; 
m1=1:r1; 
[m,m1]=meshgrid(m,m1);%生成网格空间 
noise=20.*imnoise(zeros(r,r1),'gaussian',0,0.008);%高斯噪声 
%figure(3);  
%subplot(1,2,1);
%imshow(noise,[]);
%title('白噪声') ;
a=double(21/100);%x方向的最大移动量为ra的0.21倍,可调 
b=double(21/100);%y方向的最大移动量为ca的0.21倍,可调 
t=double(88/100);%移动到最大所需的时间默认为0.88
f=ones(r,r1);   
g=(m-r/2-1).*a+(m1-r1/2-1).*b+eps;   
f=t.*sin(pi.*g).*exp(-j.*pi.*g)./(pi.*g);
h=f'.*y1; 
tu=ifft2(h);  
tu=abs(tu)+noise; 
%subplot(1,2,2);
%imshow(tu,[]);
%title('加入白噪声退化的图像')%原图傅立叶变换估计值 
axes(handles.axes2);  
imshow(tu,[]);  %显示图片
handles.img=tu;
guidata(hObject,handles);

h=fftshift(fft2(tu));
x=fftshift(fft2(noise));
K=x.*conj(x)./(y1.*conj(y1));%计算K值  
w=(f.*conj(f))'.*h./(f.*(f.*conj(f)+K'))'; 
weina=abs(ifft2(w)); 
%subplot(1,2,2);
%imshow(weina,[]);
%title('维纳滤波的结果');
axes(handles.axes3);  
imshow(weina,[]);  %显示图片
handles.img=weina;
guidata(hObject,handles);




% 有约束最小二乘滤波复原--------------------------------------------------------------------
function Untitled_4_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%I = getimage(gca);
%I=rgb2gray(I);
%I=double(I);
I = im2double(getimage(gca));
[hei,wid,~] = size(I);
%subplot(2,2,1),imshow(I);
%title('原图像');
% 模拟运动模糊.
LEN = 21;
THETA = 11;
PSF = fspecial('motion', LEN, THETA);%产生运动模糊算子，即点扩展函数
blurred = imfilter(I, PSF, 'conv', 'circular');
%subplot(2,2,2), imshow(blurred); title('模糊图像');
Pf = psf2otf(PSF,[hei,wid]);%退化函数的FFT
% 添加加性噪声
noise_mean = 0;
noise_var = 0.00001;
blurred_noisy = imnoise(blurred, 'gaussian',noise_mean, noise_var);
%subplot(2,2,3), imshow(blurred_noisy)
%title('带运动模糊和噪声图像')
axes(handles.axes2);  
imshow(blurred_noisy);  %显示图片
handles.img=blurred_noisy;
guidata(hObject,handles);

p = [0 -1 0;-1 4 -1;0 -1 0];%拉普拉斯模板
P = psf2otf(p,[hei,wid]);
gama = 0.001;
If = fft2(blurred_noisy);
numerator = conj(Pf);%conj函数，用于求一个复数的复共轭
denominator = Pf.^2 + gama*(P.^2);
deblurred2 = ifft2( numerator.*If./ denominator );%约束最小二乘方滤波在频率域中的表达式
%subplot(2,2,4), imshow(deblurred2)
%title('约束最小二乘方滤波后图像');
axes(handles.axes3);  
imshow(deblurred2);  %显示图片
handles.img=deblurred2;
guidata(hObject,handles);


% 盲卷积滤波复原--------------------------------------------------------------------
function Untitled_5_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%I = imread('peppers.tiff');
I=getimage(gca);
A=rgb2gray(I);
f = im2double(A);
%subplot(1, 3, 1), imshow(f), title('原始图像')
F = fftshift(fft2(f));

[M, N] = size(F);
[u, v] = meshgrid(1:N, 1:M);
k = 0.0025;
H = exp(-k*((v-M/2).^2+(u-N/2).^2).^(5/6));

G = F.*H;%退化图像=退化函数*原始图像
g = ifft2(ifftshift(G));
g = uint8(abs(g)*255);
%subplot(1, 3, 2), imshow(g), title('退化图像')
axes(handles.axes2); 
imshow(g);  %显示图片
handles.img=g;
guidata(hObject,handles);

I = deconv(g, H, 110);       % 可尝试不同的半径，128、108、78、48
%subplot(1, 3, 3), imshow(I), title('复原图像')
axes(handles.axes3);  
imshow(I);  %显示图片
handles.img=I;
guidata(hObject,handles);

function I_new = deconv(I, H, thresh)
if size(I, 3) == 3
    I = rgb2gray(I);
end
I = im2double(I);
G = fftshift(fft2(I));
[M, N] = size(G);
F = G;
[x, y] = meshgrid(1:N, 1:M);
if thresh > M/2
    F = G./(H+eps);
else
    idx = (x-N/2).^2 + (y-M/2).^2 < thresh^2;
    F(idx) = G(idx)./(H(idx)+eps);
end
I_new = ifft2(ifftshift(F));
I_new = uint8(abs(I_new)*255);

%Lucy-Richardson滤波复原 --------------------------------------------------------------------
function Untitled_6_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 原始图像
%f = checkerboard(8);
I=getimage(gca);
f=rgb2gray(I);
f = im2double(f);
% 噪声滤波器
PSF = fspecial('motion', 7, 45);%特殊滤波器
% 退化图像
gb = imfilter( f, PSF, 'circular' );
% 高斯滤波
noise = imnoise( zeros(size(f)), 'gaussian', 0, 0.001 );
% 将噪声加到原图上
g = gb + noise;
%figure
%subplot(221), imshow( f ), title('原图')
%subplot(222), imshow( noise, [] ), title('高斯噪声')
%subplot(223), imshow( gb ), title('退化图像')
%subplot(224), imshow( g ), title('退化图像加高斯噪声')
axes(handles.axes2);  
imshow(gb);  %显示图片
handles.img=gb;
guidata(hObject,handles);

%% Lucy-Richardson算法的迭代非线性复原
%g = checkerboard(8);
I=getimage(gca);
%f=rgb2gray(I);
g = im2double(I);
ori = g;
PSF = fspecial( 'gaussian', 7, 10 );
SD = 0.01;
g = imnoise( imfilter(g, PSF), 'gaussian', 0, SD^2 );

damper = 10*SD;
lim = ceil( size(PSF, 1) / 2 );
weight = zeros( size(g) );
weight( lim+1:end-lim, lim+1:end-lim ) = 1;
numit = 5;
f5 = deconvlucy( g, PSF, numit, damper, weight );
numit = 20;
f20 = deconvlucy( g, PSF, numit, damper, weight );
numit = 50;
f50 = deconvlucy( g, PSF, numit, damper, weight );
numit = 100;
f100 = deconvlucy( g, PSF, numit, damper, weight );
%figure
%subplot(231), imshow( ori ), title('原始图像')
%subplot(232), imshow( g, [] ), title('加两次高斯噪声的图')
%subplot(233), imshow( f5, [] ), title('LR 5次迭代')
%subplot(234), imshow( f20 ), title('LR 20次迭代')
%subplot(235), imshow( f50, [] ), title('LR 50次迭代')
%subplot(236), imshow( f100,[] ), title('LR 100次迭代')
axes(handles.axes3);  
imshow(f100, []);  %显示图片
handles.img=f100,[];%LR5次迭代
guidata(hObject,handles);

% 读入图像--- Executes on button press in pushbutton1.
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


%清除结果--- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1); %指定需要清空的坐标轴
cla reset;
axes(handles.axes2); %指定需要清空的坐标轴
cla reset;
axes(handles.axes3); %指定需要清空的坐标轴
cla reset;

%图像模糊菜单栏 --------------------------------------------------------------------
function Untitled_7_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%运动模糊 --------------------------------------------------------------------
function Untitled_8_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = im2double(getimage(gca));
[hei,wid,~] = size(I);
%subplot(2,2,1),imshow(I);
%title('原图像');
% 模拟运动模糊.
LEN = 21;
THETA = 11;
PSF = fspecial('motion', LEN, THETA);%产生运动模糊算子，即点扩展函数
blurred = imfilter(I, PSF, 'conv', 'circular');
%subplot(2,2,2), imshow(blurred); title('模糊图像');
Pf = psf2otf(PSF,[hei,wid]);%退化函数的FFT
% 添加加性噪声
noise_mean = 0;
noise_var = 0.00001;
blurred_noisy = imnoise(blurred, 'gaussian',noise_mean, noise_var);
%subplot(2,2,3), imshow(blurred_noisy)
%title('带运动模糊和噪声图像')
axes(handles.axes2);  
imshow(blurred_noisy);  %显示图片
handles.img=blurred_noisy;
guidata(hObject,handles);

% 散焦模糊--------------------------------------------------------------------
function Untitled_9_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = getimage(gca);
%I=rgb2gray(I);
I=im2double(I);
%img=imread('peppers.tiff');
%img=ind2gray(I);
%figure,imshow(img),title('原图');
psf=fspecial('disk',10);%disk掩模
res1=deconvblind(I,psf);
%figure,imshow(res1),title('盲去卷积10次');
axes(handles.axes2);  
imshow(res1);  %显示图片

%湍流模型模糊 --------------------------------------------------------------------
function Untitled_10_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% 高斯模糊--------------------------------------------------------------------
function Untitled_11_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = getimage(gca);
%imshow(I);

%% 选择噪声
m=1;    % m=1 选择高斯噪声
        % m=2 选择椒盐噪声
        
if m==1         % 添加高斯噪声
    %figure('name','添加高斯噪声','NumberTitle','off');
    J=imnoise(I,'gaussian',0,0.1);
    axes(handles.axes2);  
    imshow(J);
else if m==2    % 添加椒盐噪声
        %figure('name','添加椒盐噪声','NumberTitle','off');
        J=imnoise(I,'salt & pepper',0.1);
        %imshow(J);
        axes(handles.axes2);  
        imshow(J);  %显示图片
     end
end
