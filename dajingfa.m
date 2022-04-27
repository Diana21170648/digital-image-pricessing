%大津法
I=imread('peppers.tiff');
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
figure('name','Ostu','NumberTitle','off');
imshow(seg_I);


