p = imread('airplane.tiff');
[srcHist,srcX] = imhist(p,256);
srcHist = srcHist/numel(p);
cps = zeros(256,1,'double');
for i=1:1:256
    cps(i) = sum(srcHist(1:i));
end
cps = 255*cps;
cps = uint8(cps);

P1 = imread('peppers.tiff'); 
[dstHist,dstX] = imhist(P1,256);
dstHist = dstHist/sum(dstHist); 
cpd = zeros(256,1,'double');
for i=1:1:256
    cpd(i) = sum(dstHist(1:i));
end
cpd = cpd*255; 
cpd = uint8(cpd);

srcl=zeros(256,1,'uint8');
minv = 256;
for i = 1:1:256 
    minv =256;
    for j = 1:1:256
        if minv > abs(cps(i)-cpd(j))
            minv =  abs(cps(i)-cpd(j));
            srcl(i) =j;
        end
    end
end
[width,height] = size( p);
gray1 = p;
for i=1:1:width
    for j = 1:1:height
        gray1(i,j)=srcl(p(i,j)+1);
    end
end

[g1Hist,g1X] = imhist(gray1);
% g1Hist = g1Hist/sum(g1Hist);
g2 = histeq(p,dstHist);
[g2Hist,g2X] = imhist(g2);
% g2Hist = g2Hist/sum(g2Hist);
figure(1),
subplot(2,3,1),imshow(p),title('原始图像');
subplot(2,3,2),imshow(P1),title('匹配图像');
subplot(2,3,3),imshow(gray1),title('规定化结果图');
%subplot(2,4,4),imshow(g2),title('结果2');
subplot(2,3,4),stem(srcX,srcHist),title('原始图像直方图');
subplot(2,3,5),stem(dstX,dstHist),title('匹配图像直方图');
subplot(2,3,6),stem(g1X,g1Hist),title('规定化结果直方图');
%subplot(2,4,8),stem(g2X,g2Hist),title('结果2直方图');
