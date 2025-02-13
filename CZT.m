clear;
clc;
tic;

% 参数设置
R = 0.855;  % 目标曲率半径
lamda = 589.3e-9;  % 波长
ruler = 0.005 / 318.78;  % 计算比例尺

% 读取图像并预处理
imagefilename = "";
Ic = im2double(imread(imagefilename));
Ic = Ic - mean(Ic(:));  % 去除均值

% 设置环心坐标
m0 = 120;  % 环心x坐标
n0 = 180;  % 环心y坐标

% 图像尺寸和坐标网格
[L2, L1] = size(Ic);
[X, Y] = meshgrid(linspace(-(L1-2)/2, L1/2, L1), linspace(-(L2-2)/2, L2/2, L2));
x0 = X(1, m0);  % 圆心 x 坐标
y0 = Y(n0, 1);  % 圆心 y 坐标

% 计算每个像素点的距离
r2 = zeros(L1 * L2, 2);
for i = 1:L1
    for j = 1:L2
        r2((i-1) * L2 + j, 1) = (X(1, i) - x0)^2 + (Y(j, 1) - y0)^2;  % 像素到圆心的距离
        r2((i-1) * L2 + j, 2) = Ic(j, i);  % 对应的像素强度
    end
end

% 排序并进行采样
B_1 = sortrows(r2, 1);
C = floor(10 * numel(B_1) / 10);
B = B_1(1:C, :);

% 设置采样和频率
dt = 1 * ruler^2;
fs = 1.0 / dt;

% 信号矩阵和傅里叶变换
signal_last = reshape([B(:, 1), B(:, 2)], [], 1);
Ic_fft = fft(signal_last, numel(signal_last));

% 计算曲率半径误差
[h, g] = max(Ic_fft(1:end/2));
K_e = g * fs / numel(signal_last);
Re = 1 / (lamda * K_e);
fprintf('半径误差为 e = %6.17e\n', (abs(Re - R) / R) * 100);

% 设置 Chirp-Z 变换参数
f1 = K_e - 20000;
f2 = K_e + 20000;
m = 2048;
w = exp(-1j * 2 * pi * (f2 - f1) / fs / m);
a = exp(1j * 2 * pi * f1 / fs);

% 进行 Chirp-Z 变换
z = czt(signal_last, m, w, a);
f0 = (f2 - f1) / m * (0:m-1) + f1;

% 绘制图形
figure;
plot(f0, abs(z));
title('CZT')

% 计算并输出CZT的曲率半径误差
[aa2, bb2] = max(abs(z));
max_f2 = (f2 - f1) / m * bb2 + f1;
Re2 = 1 / (lamda * max_f2);
fprintf('CZT峰值的频率为：%6.2f Hz\n', max_f2);
fprintf('CZT半径误差为 e = %6.17e\n', (abs(Re2 - R) / R) * 100);

toc;
