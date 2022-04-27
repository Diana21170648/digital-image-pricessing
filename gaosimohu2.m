%<span style="font-size:14px;">clear
clc
%% 读入并显示原图
figure('name','原图','NumberTitle','off');
I=imread('cameraman.tif');
imshow(I);

%% 选择噪声
m=1;    % m=1 选择高斯噪声
        % m=2 选择椒盐噪声
        
if m==1         % 添加高斯噪声
    figure('name','添加高斯噪声','NumberTitle','off');
    J=imnoise(I,'gaussian',0,0.1);
    imshow(J);
else if m==2    % 添加椒盐噪声
        figure('name','添加椒盐噪声','NumberTitle','off');
        J=imnoise(I,'salt & pepper',0.1);
        imshow(J);
     end
end

%% 选择滤波器
n=2;    % n=1 选择均值滤波
        % n=2 选择中值滤波
k=1;    % 标识分隔窗口时的窗口号
w=1;    % 同时选择两种滤波器

if w==1||n==1			% MY均值滤波
    figure('name','MY均值滤波','NumberTitle','off');
    for j=3:2:9
        J1=myfilter2('average',j,J);
        subplot(2,2,k);
        subimage(J1);
        title([num2str(j),'x',num2str(j),'模版']);
        k=k+1;
    end
end
k=1;
if w==1||n==2			% MY中值滤波
    figure('name','MY中值滤波','NumberTitle','off');
    for j=3:2:9
        J2=myfilter2('medium',j,J);
        subplot(2,2,k);
        subimage(J2);
        title([num2str(j),'x',num2str(j),'模版']);
        k=k+1;
    end
end
%</span>

