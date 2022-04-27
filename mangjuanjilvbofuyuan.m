%盲卷积滤波复原
I = imread('peppers.tiff');
A=rgb2gray(I);
f = im2double(A);
subplot(1, 3, 1), imshow(f), title('原始图像')
F = fftshift(fft2(f));

[M, N] = size(F);
[u, v] = meshgrid(1:N, 1:M);
k = 0.0025;
H = exp(-k*((v-M/2).^2+(u-N/2).^2).^(5/6));

G = F.*H;
g = ifft2(ifftshift(G));
g = uint8(abs(g)*255);
subplot(1, 3, 2), imshow(g), title('退化图像')

I = deconv(g, H, 120);       % 可尝试不同的半径，128、108、78、48
subplot(1, 3, 3), imshow(I), title('复原图像')

function I_new = deconv(I, H, thresh)
if size(I, 3) == 3
    I = rgb2gray(I);
end
I = im2double(I);
G = fftshift(fft2(I));
[M, N] = size(G);
F = G;
[x, y] = meshgrid(1:N, 1:M);
if thresh > M/2,
    F = G./(H+eps);
else
    idx = (x-N/2).^2 + (y-M/2).^2 < thresh^2;
    F(idx) = G(idx)./(H(idx)+eps);
end
I_new = ifft2(ifftshift(F));
I_new = uint8(abs(I_new)*255);
end
