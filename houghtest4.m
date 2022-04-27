%%Hougn变换
I=imread('车道检测.jpg');
f=rgb2gray(I);%RGB-->gray
f=f(round(end/2):end,1:round(end/2));
BW=edge(f,'canny');  %edge:以灰度图像为输入，'canny'为边缘检测算子
                     %     输出BW为二值图像，边缘处为白（255）其余部分为黑（0）
imshow(f)
title('原始图像')
[row,col]=size(BW);
rhomax=round((row*row+col*col)^0.5);
A=zeros(2*rhomax,180);   %这里，实际上rho的取值范围为[-rhomax,rhomax]，
                         %但是为了后面进行数量统计，转变为[1,2rhomax]
for m=1:row
    for n=1:col
        if BW(m,n)>0 %判断为边缘
            for theta=1:180
                r=theta/180*pi; %角度转换
                rho=round(m*cos(r)+n*sin(r));
                %Hough变换
                rho=rho+rhomax+1;   %坐标平移
                                    %这里的理解，首先matlab中数组是从1开始计数，所以+1；
                                    %数组坐标不能<0，所以 +rhomax
                A(rho,theta)=A(rho,theta)+1;   %数量统计
            end
        end
    end
end
[rho,theta]=find(A>130);   %超过130个点视为共线,rho列号，theta行号
nma=length(rho);
figure,imshow(BW)
for i=1:nma
    hold on
    m=1:row;
    r=theta(i)/180*pi;
    n=(rho(i)-rhomax-m*cos(r))/(0.0001+sin(r));
    plot(n,m,'w-','LineWidth',6);
end
title('hough线检测');

%%MATLAB工具箱中的Hough变换函数进行边缘检测
I=imread('车道检测.jpg');
rotI=imrotate(I,33,'crop');
BW=edge(rotI,'canny');
imshow(BW);
title('原图');
%对图像进行Hough变换
[H,T,R]=hough(BW);  %[H,theta,rho]
%显示变换域
figure,imshow(imadjust(rescale(H)),'XData',T,'YData',R,...
              'InitialMagnification','fit');
xlabel('\theta');ylabel('\rho');
axis on,axis normal,hold on
title('变换域');
%计算变换域峰值
P=houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
x=T(P(:,2));y=R(P(:,1));
plot(x,y,'s','color','red');
%标记直线
lines=houghlines(BW,T,R,P,'FillGap',5,'MinLength',7);
figure,imshow(rotI),hold on
max_len=0;
for k=1:length(lines)
    xy=[lines(k).point1;lines(k).point2];
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','white');
    %Plot beginning and ends of lines
    plot(xy(1,1),xy(1,2),'xw','LineWidth',2);
    plot(xy(2,1),xy(2,2),'xw','LineWidth',2);
    %Determine the endpoints of the longest line segment
    len=norm(lines(k).point1-lines(k).point2);
    if(len>max_len)
        max_len=len;
        xy_long=xy;
    end
end
title('检测结果')