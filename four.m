function varargout = four(varargin)
% FOUR MATLAB code for four.fig
%      FOUR, by itself, creates a new FOUR or raises the existing
%      singleton*.
%
%      H = FOUR returns the handle to a new FOUR or the handle to
%      the existing singleton*.
%
%      FOUR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FOUR.M with the given input arguments.
%
%      FOUR('Property','Value',...) creates a new FOUR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before four_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to four_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help four

% Last Modified by GUIDE v2.5 19-Apr-2022 20:52:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @four_OpeningFcn, ...
                   'gui_OutputFcn',  @four_OutputFcn, ...
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

% --- Executes just before four is made visible.
function four_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to four (see VARARGIN)

% Choose default command line output for four
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using four.
if strcmp(get(hObject,'Visible'),'off')
    plot(rand(5));
end

% UIWAIT makes four wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = four_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA
% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
cla;

popup_sel_index = get(handles.popupmenu1, 'Value');
switch popup_sel_index
    case 1
        plot(rand(5));
    case 2
        plot(sin(1:0.01:25.99));
    case 3
        bar(1:.5:10);
    case 4
        plot(membrane);
    case 5
        surf(peaks);
end

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});

%理想高通 --------------------------------------------------------------------
function Untitled_8_Callback(hObject, eventdata, handles)
color_pic = getimage(gca);

gray_pic=rgb2gray(color_pic);    %将彩色图转换成灰度图
double_gray_pic=im2double(gray_pic);   %将uint8转成im2double型便于后期计算

[width,height]=size(double_gray_pic);
mid_w=width/2;    %图像中心点横坐标
mid_h=height/2;   %图像中心点纵坐标
fourier_pic=fft2(double_gray_pic);   %对灰度图进行傅里叶变换
g=fftshift(fourier_pic);  %将频谱图中零频率成分移动至频谱图中心

end_radius=[10,30,80,150];   %设置截止频率
%Result=zeros(width,height);  %预先分配内存空间，提高运行速率
%figure('name','理想低通滤波器');
%subplot(3,2,1);imshow(double_gray_pic,[]);title('原灰度图');
for i=1:width
        for j=1:height
            distance=sqrt((i-mid_w)^2+(j-mid_h)^2);   %计算点（x,y）到中心点的距离
              if distance>end_radius(1)  %如果距离大于截止频率，则滤除分量，直接置0
                Result1(i,j)=g(i,j);
              end
        end
    end
      output1=im2uint8(real(ifft2(ifftshift(Result1))));  
      
        for i=1:width
        for j=1:height
            distance=sqrt((i-mid_w)^2+(j-mid_h)^2);   %计算点（x,y）到中心点的距离
              if distance>end_radius(2)  %如果距离大于截止频率，则滤除分量，直接置0
                Result2(i,j)=g(i,j);%这样会陷入死循环，不弄%result2(i,j)=fourier_shift(i,j)*h2;
              end
        end
    end
   output2=im2uint8(real(ifft2(ifftshift(Result2)))); 
    %subplot(3,2,k+2);imshow(output,[]);title(['理想低通滤波器 D0=',num2str(end_radius(k))]); 

for i=1:width
        for j=1:height
            distance=sqrt((i-mid_w)^2+(j-mid_h)^2);   %计算点（x,y）到中心点的距离
              if distance>end_radius(3)  %如果距离大于截止频率，则滤除分量，直接置0
                Result3(i,j)=g(i,j);
              end
        end
    end
      output3=im2uint8(real(ifft2(ifftshift(Result3))));  
      
        for i=1:width
        for j=1:height
            distance=sqrt((i-mid_w)^2+(j-mid_h)^2);   %计算点（x,y）到中心点的距离
              if distance>end_radius(4)  %如果距离大于截止频率，则滤除分量，直接置0
                Result4(i,j)=g(i,j);%这样会陷入死循环，不弄%result2(i,j)=fourier_shift(i,j)*h2;
              end
        end
    end
      output4=im2uint8(real(ifft2(ifftshift(Result4))));  
    %subplot(3,2,k+2);imshow(output,[]);title(['理想低通滤波器 D0=',num2str(end_radius(k))]); 

axes(handles.axes2);
imshow(output1);
handles.img=output1;
guidata(hObject,handles);

axes(handles.axes3);
imshow(output2);
handles.img=output2;
guidata(hObject,handles);

axes(handles.axes4);
imshow(output3);
handles.img=output3;
guidata(hObject,handles);

axes(handles.axes5);
imshow(output4);
handles.img=output4;
guidata(hObject,handles);

% 巴特沃斯高通滤波器--------------------------------------------------------------------
function Unt_9_Callback(hObject, eventdata, handles)
color_pic= getimage(gca);

gray_pic=rgb2gray(color_pic);    %将彩色图转换成灰度图
double_gray_pic=im2double(gray_pic);   %将uint8转成im2double型便于后期计算

[width,height]=size(double_gray_pic);
mid_w=width/2;    %图像中心点横坐标
mid_h=height/2;   %图像中心点纵坐标
fourier_pic=fft2(double_gray_pic);   %对灰度图进行傅里叶变换
fourier_shift=fftshift(fourier_pic);  %将频谱图中零频率成分移动至频谱图中心

level=2;   %二阶巴特沃兹
end_radius=[10,30,80,150];    %设置截止频率
result1=zeros(width,height);  %预先分配内存空间，提高运行速率
result2=zeros(width,height);  %预先分配内存空间，提高运行速率
result3=zeros(width,height);  %预先分配内存空间，提高运行速率
for i=1:width
    for j=1:height
        distance=sqrt((i-mid_w)^2+(j-mid_h)^2);   %计算点（x,y）到中心点的距离
        h1=1./(1+(end_radius(1)/distance).^(2*level)); %计算巴特沃斯滤波器,H = mat2gray(1./(1+((D./D0).^4)));
        h2=1./(1+(end_radius(2)/distance).^(2*level));
        h3=1./(1+(end_radius(3)/distance).^(2*level));
               h4=1./(1+(end_radius(4)/distance).^(2*level));
 
        result1(i,j)=fourier_shift(i,j)*h1;   %用滤波器乘以主函数
        result2(i,j)=fourier_shift(i,j)*h2;
        result3(i,j)=fourier_shift(i,j)*h3;
                result4(i,j)=fourier_shift(i,j)*h4;

    end
end 
output1=im2uint8(real(ifft2(ifftshift(result1))));  %最终输出要记得频谱搬移回去
output2=im2uint8(real(ifft2(ifftshift(result2))));
output3=im2uint8(real(ifft2(ifftshift(result3))));
output4=im2uint8(real(ifft2(ifftshift(result4))));

axes(handles.axes2);
imshow(output1);
handles.img=output1;
guidata(hObject,handles);

axes(handles.axes3);
imshow(output2);
handles.img=output2;
guidata(hObject,handles);

axes(handles.axes4);
imshow(output3);
handles.img=output3;
guidata(hObject,handles);

axes(handles.axes5);
imshow(output4);
handles.img=output4;
guidata(hObject,handles);

%figure('name','巴特沃兹高通滤波器');
%subplot(3,2,1);imshow(double_gray_pic);title('原灰度图');
%subplot(3,2,3);imshow(output1,[]);title(['巴特沃兹高通滤波 D0=',num2str(end_radius(1))]);
%subplot(3,2,4);imshow(output2,[]);title(['巴特沃兹高通滤波 D0=',num2str(end_radius(2))]);
%subplot(3,2,5);imshow(output3,[]);title(['巴特沃兹高通滤波 D0=',num2str(end_radius(3))]);
%subplot(3,2,6);imshow(output4,[]);title(['巴特沃兹高通滤波 D0=',num2str(end_radius(4))]);

% 高斯高通--------------------------------------------------------------------
function Untitled_10_Callback(hObject, eventdata, handles)
imgrgb = getimage(gca);

f = rgb2gray(imgrgb); %将rgb图像转换成灰度图像
 
%高斯高通滤波
I = double(f);
g = fft2(I);%二维傅立叶变换
g = fftshift(g);%频移
[M, N] = size(g);
D0 = [10,30,80,150];%截止频率为5
m = fix(M / 2); n = fix(N / 2);
 
for i = 1:M
    for j = 1:N
        D = sqrt((i - m)^2 + (j - n)^2);
        H1 = exp(-(D.^2) ./ (2 * (D0(1)^2)));
        result1(i, j) = (1- H1) * g(i, j);
    end
end
result1 = ifftshift(result1);
J1 = ifft2(result1);
J2 = uint8(real(J1));
%subplot(3, 2, 3);
%imshow(J2)
%title('高斯高通滤波 D0=10')
axes(handles.axes2);
imshow(J2);%(output4);
handles.img=J2;%output4;
guidata(hObject,handles)

for i = 1:M
 
    for j = 1:N
        D = sqrt((i - m)^2 + (j - n)^2);
        H2 = (exp(-(D.^2) ./ (2 * (D0(2)^2))));
        result2(i, j) =  (1- H2) * g(i, j);
    end
end
 
result2 = ifftshift(result2);
J1 = ifft2(result2);
J2 = uint8(real(J1));
%subplot(3, 2, 4);
%imshow(J2)
%title('高斯高通滤波 D0=30')
axes(handles.axes3);
imshow(J2);%(output4);
handles.img=J2;%output4;
guidata(hObject,handles)

for i = 1:M
 
    for j = 1:N
        D = sqrt((i - m)^2 + (j - n)^2);
        H3 = exp(-(D.^2) ./ (2 * (D0(3)^2)));
        result3(i, j) =  (1- H3) * g(i, j);
    end
end 
result3 = ifftshift(result3);
J1 = ifft2(result3);
J2 = uint8(real(J1));
%subplot(3, 2, 5);
%imshow(J2)
%title('高斯高通滤波 D0=80')
axes(handles.axes4);
imshow(J2);%(output4);
handles.img=J2;%output4;
guidata(hObject,handles)

for i = 1:M
 
    for j = 1:N
        D = sqrt((i - m)^2 + (j - n)^2);
        H4 = exp(-(D.^2) ./ (2 * (D0(4)^2)));
        result4(i, j) =(1- H4) * g(i, j);
    end
end
 
result4 = ifftshift(result4);
J1 = ifft2(result4);
J2 = uint8(real(J1));
%subplot(3, 2, 6);
%imshow(J2)
%title('高斯高通滤波 D0=150')
axes(handles.axes5);
imshow(J2);%(output4);
handles.img=J2;%output4;
guidata(hObject,handles)

%理想低通滤波
function U_3_Callback(hObject, eventdata, handles)
%理想低通
I = getimage(gca);
I=rgb2gray(I);
%figure(1);
%subplot(321),imshow(I);
%title('原图像');
%I=imnoise(I,'gaussian');%%加入高斯白噪声
%subplot(222),imshow(I);
%title('加入噪声后的图像');
s=fftshift(fft2(I));
%subplot(223), imshow(log(abs(s)),[]); 
%title('图像傅里叶变换取对数所得频谱');
[a,b]=size(s);
a0=round(a/2);
b0=round(b/2);
d=[10,30,80,150];
for i=1:a 
    for j=1:b 
        distance=sqrt((i-a0)^2+(j-b0)^2);
        if distance<=d(1)
            h=1;
        else
            h=0;
        end
        s1(i,j)=h*s(i,j);
    end
end
output1=uint8(real(ifft2(ifftshift(s1))));
%subplot(323),imshow(output1);
%title('10理想低通滤波所得图像'); 

for i=1:a 
    for j=1:b 
        distance=sqrt((i-a0)^2+(j-b0)^2);
        if distance<=d(2)
            h=1;
        else
            h=0;
        end
        s2(i,j)=h*s(i,j);
    end
end
output2=uint8(real(ifft2(ifftshift(s2))));
%subplot(324),imshow(output2);
%title('30理想低通滤波所得图像'); 

for i=1:a 
    for j=1:b 
        distance=sqrt((i-a0)^2+(j-b0)^2);
        if distance<=d(3)
            h=1;
        else
            h=0;
        end
        s3(i,j)=h*s(i,j);
    end
end
output3=uint8(real(ifft2(ifftshift(s3))));
%subplot(325),imshow(output3);
%title('80理想低通滤波所得图像'); 

for i=1:a 
    for j=1:b 
        distance=sqrt((i-a0)^2+(j-b0)^2);
        if distance<=d(4)
            h=1;
        else
            h=0;
        end
        s4(i,j)=h*s(i,j);
    end
end
output4=uint8(real(ifft2(ifftshift(s4))));
%subplot(326),imshow(output4);
%title('150理想低通滤波所得图像'); 


axes(handles.axes2);
imshow(output1);
handles.img=output1;
guidata(hObject,handles);

axes(handles.axes3);
imshow(output2);
handles.img=output2;
guidata(hObject,handles);

axes(handles.axes4);
imshow(output3);
handles.img=output3;
guidata(hObject,handles);

axes(handles.axes5);
imshow(output4);
handles.img=output4;
guidata(hObject,handles);

% 巴特沃斯低通--------------------------------------------------------------------
function Untitled_5_Callback(hObject, eventdata, handles)
color_pic = getimage(gca);
gray_pic=rgb2gray(color_pic);    %将彩色图转换成灰度图
double_gray_pic=im2double(gray_pic);   %将uint8转成im2double型便于后期计算
[width,height]=size(double_gray_pic);
mid_w=width/2;    %图像中心点横坐标
mid_h=height/2;   %图像中心点纵坐标
fourier_pic=fft2(double_gray_pic);   %对灰度图进行傅里叶变换
fourier_shift=fftshift(fourier_pic);  %将频谱图中零频率成分移动至频谱图中心

level=2;   %二阶巴特沃兹
end_radius=[10,30,80,150];    %设置截止频率
result1=zeros(width,height);  %预先分配内存空间，提高运行速率
result2=zeros(width,height);  %预先分配内存空间，提高运行速率
result3=zeros(width,height);  %预先分配内存空间，提高运行速率
result4=zeros(width,height);  %预先分配内存空间，提高运行速率
for i=1:width
    for j=1:height
        distance=sqrt((i-mid_w)^2+(j-mid_h)^2);   %计算点（x,y）到中心点的距离
        h1=1./(1+(distance/end_radius(1)).^(2*level)); %计算巴特沃斯滤波器
        h2=1./(1+(distance/end_radius(2)).^(2*level));
        h3=1./(1+(distance/end_radius(3)).^(2*level));
        h4=1./(1+(distance/end_radius(4)).^(2*level));
        result1(i,j)=fourier_shift(i,j)*h1;   %用滤波器乘以主函数,H = mat2gray(1./(1+((D./D0).^4)));
        result2(i,j)=fourier_shift(i,j)*h2;
        result3(i,j)=fourier_shift(i,j)*h3;
        result4(i,j)=fourier_shift(i,j)*h4;
    end
end
output1=im2uint8(real(ifft2(ifftshift(result1))));  %最终输出要记得频谱搬移回去
output2=im2uint8(real(ifft2(ifftshift(result2))));
output3=im2uint8(real(ifft2(ifftshift(result3))));
output4=im2uint8(real(ifft2(ifftshift(result4))));

axes(handles.axes2);
imshow(output1);
handles.img=output1;
guidata(hObject,handles);

axes(handles.axes3);
imshow(output2);
handles.img=output2;
guidata(hObject,handles);

axes(handles.axes4);
imshow(output3);
handles.img=output3;
guidata(hObject,handles);

axes(handles.axes5);
imshow(output4);
handles.img=output4;
guidata(hObject,handles);

%figure('name','巴特沃兹低通滤波器');
%subplot(3,2,1);imshow(double_gray_pic);title('原灰度图');
%subplot(3,2,3);imshow(output1,[]);title(['巴特沃兹低通滤波 D0=',num2str(end_radius(1))]);
%subplot(3,2,4);imshow(output2,[]);title(['巴特沃兹低通滤波 D0=',num2str(end_radius(2))]);
%subplot(3,2,5);imshow(output3,[]);title(['巴特沃兹低通滤波 D0=',num2str(end_radius(3))]);
%subplot(3,2,6);imshow(output4,[]);title(['巴特沃兹低通滤波 D0=',num2str(end_radius(4))]);

%高斯低通 --------------------------------------------------------------------
function Untitled_6_Callback(hObject, eventdata, handles)
imgrgb = getimage(gca);

f = rgb2gray(imgrgb); %将rgb图像转换成灰度图像
 
%高斯低通滤波
I = double(f);
g = fft2(I);%二维傅立叶变换
g = fftshift(g);%频移
[M, N] = size(g);
D0 = [10,30,80,150];%截止频率为5
m = fix(M / 2); n = fix(N / 2);
 
for i = 1:M
    for j = 1:N
        D = sqrt((i - m)^2 + (j - n)^2);
        H1 = exp(-(D.^2) ./ (2 * (D0(1)^2)));
        result1(i, j) =  H1 * g(i, j);
    end
end
result1 = ifftshift(result1);
J1 = ifft2(result1);
J2 = uint8(real(J1));
%subplot(3, 2, 3);
%imshow(J2)
%title('高斯低通滤波 D0=10')
axes(handles.axes2);
imshow(J2);%(output4);
handles.img=J2;%output4;
guidata(hObject,handles)

for i = 1:M
 
    for j = 1:N
        D = sqrt((i - m)^2 + (j - n)^2);
        H2 = (exp(-(D.^2) ./ (2 * (D0(2)^2))));
        result2(i, j) =  H2 * g(i, j);
    end
end
 
result2 = ifftshift(result2);
J1 = ifft2(result2);
J2 = uint8(real(J1));
%subplot(3, 2, 4);
%imshow(J2)
%title('高斯低通滤波 D0=30')
axes(handles.axes3);
imshow(J2);%(output4);
handles.img=J2;%output4;
guidata(hObject,handles)

for i = 1:M
 
    for j = 1:N
        D = sqrt((i - m)^2 + (j - n)^2);
        H3 = exp(-(D.^2) ./ (2 * (D0(3)^2)));
        result3(i, j) =  H3 * g(i, j);
    end
end 
result3 = ifftshift(result3);
J1 = ifft2(result3);
J2 = uint8(real(J1));
%subplot(3, 2, 5);
%imshow(J2)
%title('高斯低通滤波 D0=80')
axes(handles.axes4);
imshow(J2);%(output4);
handles.img=J2;%output4;
guidata(hObject,handles)

for i = 1:M
 
    for j = 1:N
        D = sqrt((i - m)^2 + (j - n)^2);
        H4 = exp(-(D.^2) ./ (2 * (D0(4)^2)));
        result4(i, j) = H4 * g(i, j);
    end
end
 
result4 = ifftshift(result4);
J1 = ifft2(result4);
J2 = uint8(real(J1));
%subplot(3, 2, 6);
%imshow(J2)
%title('高斯低通滤波 D0=150')
axes(handles.axes5);
imshow(J2);%(output4);
handles.img=J2;%output4;
guidata(hObject,handles)



% 傅里叶变换--------------------------------------------------------------------
function Untitled_16_Callback(hObject, eventdata, handles)
%GUI2; % 第二个界面（主界面）————推荐直接输入函数名的这种方式% 或者：
%roberts

% --------------------------------------------------------------------
function Untitled_17_Callback(hObject, eventdata, handles)
structure with handles and user data (see GUIDATA)
GUI2;

%系统读入图像
function pushbutton2_Callback(hObject, eventdata, handles)
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
global T
T=im1
axes(handles.axes1);  
imshow(im1);  %显示图片
handles.img=im1;
guidata(hObject,handles);

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
axes(handles.axes1); %指定需要清空的坐标轴
cla reset;
axes(handles.axes2); %指定需要清空的坐标轴
cla reset;
axes(handles.axes3); %指定需要清空的坐标轴
cla reset;
axes(handles.axes4); %指定需要清空的坐标轴
cla reset;
axes(handles.axes5); %指定需要清空的坐标轴
cla reset;


% 低通滤波按钮不能删--------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% 高通滤波按钮不能删--------------------------------------------------------------------
function Untitled_7_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% 图像滤波菜单不能删--------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%傅里叶变换菜单不能删--------------------------------------------------------------------
function Untitled_13_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%图像分割边缘检测打开GUI3 --------------------------------------------------------------------
function Untitled_20_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
open GUI3.fig


%复原按钮图像 --------------------------------------------------------------------
function Untitled_21_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
open GUI_4.fig
