clear;
close all;
clc;
tic;

% 参数设置
imagefilename = "";
Ic = im2double(rgb2gray(imread(imagefilename)));
Ic = Ic - mean(Ic(:));  % 去除均值

% 设置分块参数
num = 90;
[L2, L1] = size(Ic);
sub_wd = floor(L1 / num);
sub_hd = floor(L2 / num);

% 生成坐标网格
[X, Y] = meshgrid(1:L1, 1:L2);

% 窗函数
win = hann(num) * hann(num)';

% 分块傅里叶变换
k = 0;
for ii = 1:sub_hd
    for jj = 1:sub_wd
        temp = Ic((ii-1)*sub_hd+1:ii*sub_hd, (jj-1)*sub_wd+1:jj*sub_wd);
        F_C = fft2(temp .* win);
        FS_C = fftshift(F_C);  % 将低频转移到中心
        S_C = log(abs(FS_C));  % 对数变换显示傅里叶谱
        S_C(1,:) = 0;  % 去除高频
        S_C(:,1) = 0;

        % 计算傅里叶变换中心坐标
        x_0 = (jj-1)*sub_wd + sub_wd / 2;
        y_0 = (ii-1)*sub_hd + sub_hd / 2;

        % 计算方向
        [row_num, col_num] = size(S_C);
        [res_x, res_y] = find(S_C > 1);  % 找到非零点
        avg_x = mean(res_x);
        avg_y = mean(res_y);
        
        % 计算惯性矩阵的特征值，得到方向角
        matrix = [res_x - avg_x, res_y - avg_y];
        matrix2 = matrix' * matrix;
        [EigMatrix, EigValue] = eig(matrix2);
        k_a = EigMatrix(1,1);
        k_b = EigMatrix(2,1);

        theta(ii,jj) = -atan(-k_b/k_a);
    end
end

% 计算牛顿环的中心
S11 = sin(theta).^2;
S22 = cos(theta).*sin(theta);
S33 = x_0 .* S11 - y_0 .* S22;
S44 = cos(theta).^2;
S55 = x_0 .* S22 - y_0 .* S44;
S1 = sum(S11(:));
S2 = sum(S22(:));
S3 = sum(S33(:));
S4 = sum(S44(:));
S5 = sum(S55(:));

x_c = (S2 * S5 - S3 * S4) / (S2^2 - S1 * S4);
y_c = (S1 * S5 - S2 * S3) / (S2^2 - S1 * S4);

x_c = round(x_c);
y_c = round(y_c);

fprintf('圆心坐标为 [%d, %d]\n', x_c, y_c);

toc;
