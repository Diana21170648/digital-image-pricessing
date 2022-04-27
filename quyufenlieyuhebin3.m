%区域分裂与合并
clear all

I=imread('cameraman.tif');

S=qtdecomp(I,.2);

blocks=repmat(uint8(0),size(S));

for dim=[512 256 128 64 32 16 8 4 2 1];

numblocks=length(find(S==dim));

if(numblocks>0)

values=repmat(uint8(1),[dim dim numblocks]);

values(2:dim,2:dim,:)=0;

blocks=qtsetblk(blocks,S,dim,values);

end

end

figure(1);

subplot(1,2,1);

imshow(I);

title('原始图像');

subplot(1,2,2);

imshow(blocks,[]);

title('分解图像');
