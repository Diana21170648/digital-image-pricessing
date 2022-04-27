I=imread('peppers.tiff');
I=double(I);
minValue=min(min(I));
maxValue=max(max(I));
[row,col]=size(I);
Th=minValue+1;
perfactValue=10000000000;
for m=minValue+1:maxValue-1  %m=minValue+1:maxValue-1  %
    k1=1;
    k2=1;
    for i=1:row
        for j=1:col 
            if I(i,j)<m
                C1(1,k1)=I(i,j);%C1类
                k1=k1+1;
            else
                C2(1,k2)=I(i,j);%C2类
                k2=k2+1;
           
            end
        end
    end
  %求每一类的均值，方差，概率分布
    average1=mean(C1);
    variance1=0;
    for i=1:k1-1 
        variance1=variance1+(C1(1,i)-average1)^2;
    end
    variance1=variance1/(k1-1);
    p1=(k1-1)/(row*col);%分布概率,row*col是像素个数
  %求每一类的均值，方差，概率分布    
    average2=mean(C2);
    variance2=0;
    for i=1:k2-1 
        variance2=variance2+(C2(1,i)-average2)^2;%方差
    end
    p2=(k2-1)/(row*col);%分布概率
    variance2=variance2/(k2-1);
    
    newValue=p1*variance1+p2*variance2;
    if (newValue<perfactValue)
        Th=m;
        perfactValue=newValue;
    end
end
for i=1:row
    for j=1:col
        if I(i,j)>=Th
            G(i,j)=255;
        else
            G(i,j)=0;
        end
    end
end
subplot(122),imshow(G);title=('segmentation')


