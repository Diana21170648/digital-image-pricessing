%理想高通 --------------------------------------------------------------------
function Untitled_8_Callback(hObject, eventdata, handles)
color_pic=imread('airplane.tiff');  %读取彩色图像
gray_pic=rgb2gray(color_pic);    %将彩色图转换成灰度图
double_gray_pic=im2double(gray_pic);   %将uint8转成im2double型便于后期计算
[width,height]=size(double_gray_pic);
mid_w=width/2;    %图像中心点横坐标
mid_h=height/2;   %图像中心点纵坐标
fourier_pic=fft2(double_gray_pic);   %对灰度图进行傅里叶变换
fourier_shift=fftshift(fourier_pic);  %将频谱图中零频率成分移动至频谱图中心

end_radius=[10,30,80,150];   %设置截止频率
Result=zeros(width,height);  %预先分配内存空间，提高运行速率
subplot(3,2,1);imshow(double_gray_pic,[]);title('原灰度图');
for k=1:4
    Result=fourier_shift;
    for i=1:width
        for j=1:height
            distance=sqrt((i-mid_w)^2+(j-mid_h)^2);   %计算点（x,y）到中心点的距离
            if distance<end_radius(k)  %如果距离大于截止频率，则滤除分量，直接置0
                Result(i,j)=0;
            end
        end
    end
    output1=im2uint8(real(ifft2(ifftshift(Result))));  %最终输出要记得频谱搬移回去
subplot(3,2,k+2);imshow(output,[]);title(['理想高通滤波器 D0=',num2str(end_radius(k))]); 
end
