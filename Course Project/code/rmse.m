clc; clear all; close all;
ratio = [0.2, 0.25, 0.3, 0.4, 0.5, 0.6]; % fration of measurements taken from 8x8 patches
x = imread('../data/peppers.png'); % test image
x = imresize(x, 0.25); % resizing it to increase computational speed
x = double(rgb2gray(x));
rmsee1 = []; % stores RMSE for different sensing matrices
x11 = load('../mat_files/recons_peppers_1.mat');
x21 = load('../mat_files/recons_peppers_2.mat');
x31 = load('../mat_files/recons_peppers_3.mat');
x41 = load('../mat_files/recons_peppers_4.mat');
x51 = load('../mat_files/recons_peppers_5.mat');
x61 = load('../mat_files/recons_peppers_6.mat');
x_est1 = {x11, x21, x31, x41, x51, x61};
rmsee2 = []; % stores RMSE for different sensing matrices
x12 = load('../mat_files/recons_peppers_cs_1.mat');
x22 = load('../mat_files/recons_peppers_cs_2.mat');
x32 = load('../mat_files/recons_peppers_cs_3.mat');
x42 = load('../mat_files/recons_peppers_cs_4.mat');
x52 = load('../mat_files/recons_peppers_cs_5.mat');
x62 = load('../mat_files/recons_peppers_cs_6.mat');
x_est2 = {x12, x22, x32, x42, x52, x62};
for i=1:size(x_est1, 2)
    xi1 = cell2mat(x_est1(i));
    xi2 = cell2mat(x_est2(i));
    rmsee1 = [rmsee1, norm(xi1.x_est(:)-x(:))/norm(x(:))];
    rmsee2 = [rmsee2, norm(xi2.x_est(:)-x(:))/norm(x(:))];
end
figure();
plot(ratio, rmsee1, 'bo-');
hold on;
plot(ratio, rmsee2, 'ro-'); 
xlabel('ratio');
ylabel('RMSE');
title('SCS vs CS for Peppers');
legend(["SCS", "CS"]);