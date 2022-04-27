    % image_edge.m
    % 探究对真彩色图像的边缘检测方法
    % ----两种边缘检测方法(均以Laplace算子为例)
    % -----------------1. 真彩色图像->灰度图->边缘检测
    % -----------------2. 真彩色图像->分解为RGB分量->对分量进行边缘转换->对转换结果进行合并->剔除不合适的分量

    %初始化运行环境
    close all; clear; clc;
    warning off;

    % 读取图像并适当进行压缩
    I_origin = imread('airplane.tiff');

    [size_x, size_y, size_z] = size(I_origin);
    if size_x > 1080
        I_origin2 = imresize(I_origin, 1080 / double(size_x));
    else 
        I_origin2 = I_origin;
    end

   % figure(1), imshow(I_origin2), title('原图');

    % ---- 方法一
    I_gray = rgb2gray(I_origin2);
    %figure('Name', '对灰度图的边缘检测');
    %subplot(1, 2, 1), imshow(I_origin2), title('原图');
    %subplot(1, 2, 2), imshow(I_gray), title('灰度图');

    Edge_gray = edge(I_gray, 'log');

    % ---- 方法二
    % ---------------提取RGB分量并显示
    Instance_R = I_origin2(:, :, 1);
    Instance_G = I_origin2(:, :, 2);
    Instance_B = I_origin2(:, :, 3);
    figure('Name', '原图的RGB分量');
    subplot(2, 2, 1), imshow(I_origin2), title('Origin');
    subplot(2, 2, 2), imshow(Instance_R), title('Vector R');
    subplot(2, 2, 3), imshow(Instance_G), title('Vector G');
    subplot(2, 2, 4), imshow(Instance_B), title('Vector B');

    % ---------------对RGB分量进行边缘检测并合并
    Edge_R = edge(Instance_R, 'log');
    Edge_G = edge(Instance_G, 'log');
    Edge_B = edge(Instance_B, 'log');
    rgb = im2uint8(cat(3, Edge_R, Edge_G, Edge_B));

    figure('Name', 'RGB分量的边缘检测');
    subplot(2, 2, 1), imshow(I_origin2), title('Origin');
    subplot(2, 2, 2), imshow(Edge_R), title('Laplace Vector R');
    subplot(2, 2, 3), imshow(Edge_G), title('Laplace Vector G');
    subplot(2, 2, 4), imshow(Edge_B), title('Laplace Vector B');

   % figure('Name', '两种检测方法的对比');
   % subplot(1, 2, 1), imshow(Edge_gray), title('方法一');
   % subplot(1, 2, 2), imshow(rgb), title('方法二');

    %EOF
