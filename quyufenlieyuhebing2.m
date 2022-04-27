%% 使用区域分离和合并的图像分割
f = imread('peppers.tiff');
 figure,imshow(f);
 title('区域分割原始图像');
g=splitmerge(f,2,@predicate);%32代表分割中允许最小的块，predicate函数返回1，说明需要再分裂，返回0说明不需要继续分裂
figure,imshow(g);
title('mindim为2时的分割图像');
se=ones(8,8);
gdilate=imdilate(g,se);%膨胀是为了填充空洞
figure;imshow(gdilate);
title('膨胀后的图')
gerode=imerode(gdilate,se);%腐蚀是为了缩回原来大小
 figure;imshow(gerode);
 title('腐蚀后的图')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g=splitmerge(f,mindim,fun)%f是待分割的原图，mindim是定义分解中所允许的最小的块，必须是2的正整数次幂
 Q=2^nextpow2(max(size(f)));
 [M,N]=size(f);
 f=padarray(f,[Q-M,Q-N],'post');
 S=qtdecomp(f,@split_test,mindim,fun);

 Lmax=full(max(S(:)));%将以稀疏矩阵存储形式存储的矩阵变换成以普通矩阵（full matrix）形式存储，full，sparse只是存储形式的不同
 g=zeros(size(f));
 MARKER=zeros(size(f));
 
 for k=1:Lmax
     [vals,r,c]=qtgetblk(f,S,k);%vals是一个数组，包含f的四叉树分解中大小为k*k的块的值，是一个k*k*个数的矩阵，
     %个数是指S中有多少个这样大小的块，f是被四叉树分的原图像，r，c是对应的左上角块的坐标如2*2块，代表的是左上角开始块的坐标
     
     if ~isempty(vals)
         for I=1:length(r)
             xlow=r(I);
             ylow=c(I);
             xhigh=xlow+k-1;
             yhigh=ylow+k-1;
             region=f(xlow:xhigh,ylow:yhigh);%找到对应的区域
             flag=feval(fun,region);%evaluates the function handle, fhandle,using arguments x1 through xn.执行函数fun，region是参数
             if flag%如果返回的是1，则进行标记
                 g(xlow:xhigh,ylow:yhigh)=1;%然后将对应的区域置1
                 MARKER(xlow,ylow)=1;%MARKER矩阵对应的左上角坐标置1
             end
         end
     end
 end
 
 g=bwlabel(imreconstruct(MARKER,g));%imreconstruct默认2D图像8连通，这个函数就是起合的作用
 g=g(1:M,1:N);%返回原图像的大小
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function v=split_test(B,mindim,fun)
 K=size(B,3);
 v(1:K)=false;
 for I=1:K
     quadregion=B(:,:,I);
     if size(quadregion,1)<=mindim%如果分的块的大小小于mindim就直接结束
         v(I)=false;
         continue
     end
     flag=feval(fun,quadregion);%quadregion是fun函数的参数
     if flag%如果flag是1，代表需要再分
         v(I)=true;%这里就相当于split_test是起一个调用predicate的作用，返回的就是ppredicate的值
     end
 end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function flag=predicate(region)
sd=std2(region);
m=mean2(region);
flag=(sd>20)&(m>26)&(m<255);
end