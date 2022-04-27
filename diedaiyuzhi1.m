I=imread('peppers.tiff');
%I=rgb2gray(I);
I=im2double(I);
level=0.1;
max_I=max(I(:));
min_I=min(I(:));
T1=1/2*(max_I+min_I);
[M,N]=size(I);
A_number=0;
B_number=0;
A_all=0;
B_all=0;
for i=1:M
    for j=1:N
        if(I(i,j)>=T1)
            A_number=A_number+1;
            A_all=A_all+I(i,j);
        else if(I(i,j)<T1)
             B_number=B_number+1;
             B_all=B_all+I(i,j);
            end
        end
    end
    A_ave= A_all/A_number;
    B_ave= B_all/B_number;
    T2=1/2*( A_ave+ B_ave );%*不能省略
end
T2*255
J=im2bw(I,T2);
figure(1);
subplot(121);imshow(I);
subplot(122);imshow(J);
                
                
            
            