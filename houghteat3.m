Img  = imread('车道检测.jpg');
Img_gray  = rgb2gray(Img);
Img_sobel  = mysobel(Img_gray);
function img = mysobel(img)

Sobel1 = [ 1  2  1;
           0  0  0;
          -1 -2 -1];
      
Sobel2 = [-1  0  1;
          -2  0  2;
          -1  0  1];

img1 = conv2(img,Sobel1,'same');
img2 = conv2(img,Sobel2,'same');
img = sqrt(img1.^2+img2.^2);
img = (img - min(img(:)))./(max(img(:))-min(img(:)));
Img_bin  = Binpro(Img_sobel);%图像二值化
function Bin = Binpro(img)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  调参rate，通过设置rate
%%%  可以很好的去除不必要的信息
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rate = 0.98;

rank_list = sort(img(:));

listlen = length(rank_list);

t = rank_list(fix(rate*listlen));

[w ,h] = size(img);

Bin = img(i,j) >= t;

r = x*cos('\theta')+y*sin('\theta')
y = -x*cos('\theta')/sin('\theta')+r0/sin('\theta')
%建立Hough空间将点对应曲线进行映射
f = 361;

%初始化Hough空间
Hough = zeros(f,fix(sqrt(m^2+n^2))+1); 

for i = 1:m
    for j = 1:n
        if Img_bin(i,j) ~= 0
            for a = 1:f
            	
            	%%%  对 r = x*cos(θ)+y*sin(θ) 公式的使用
                a0  = a*2*pi/f;
                r = i*cos(a0) + j*sin(a0);
                
                if r >0 
                	%%%  在0-360的角度变化下，极径r的值应该是大于0的
                    r = fix(r)+1;
                    Hough(a,r) = Hough(a,r) + 1; %%通过遍历绘制曲线
                    
                end
            end
        end
    end
end

Hough0 = Hough;
Hough = uint8(Hough); %%% 为了显示Hough空间曲线绘制好的图，我们将它转化为可显示的类型，仅仅是为了显示

%5
[M,N] = size(Hough);
H.r = [];
H.a = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  调参cut_point
%%%  这个是选取在Hough空间高幅值点的标准
%%%  这里选取的是幅值从小到大排序，后0.1%的一些点点
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cut_point = 0.9999;
temp = sort(Hough0(:));
cut = temp(fix(length(temp)*cut_point));
%取出符合要求的点的r，a
for i = 1:M
    for j = 1:N
        if Hough0(i,j)>cut && compwith(i,j,H.a,H.r,20000)
            H.a = [H.a,i];
            H.r = [H.r,j];
        end
    end
end

function compwith = compwith(x,y,listx,listy,focus)

temp = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  调参focus
%%%  这是决定邻域半径的参数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(listx)
    if (listx(i) - x)^2 + (listy(i) - y)^2 < focus
        temp = temp +1;
    end
end
if temp >0
    compwith = 0;
else
    compwith = 1;
end


for i = 1:length(H.a)
    for i0 =1:m
    	a0 = H.a(i)*pi/180;
        j0 = -i0*cos(a0)/sin(a0)+H.r(i)/sin(a0);
        j0 = ceil(j0)+1;
        if j0 < n && j0 > 0  %% 约束 y的取值在矩阵的索引内
            for z = 1:4
            	%绘制显眼的线，这里用红色的
                Img(i0+z,j0,1)= 255;
                Img(i0+z,j0,2)= 0;
                Img(i0+z,j0,3)= 0;
            end
        end
    end
end
imshow(Img)
end
end
end


