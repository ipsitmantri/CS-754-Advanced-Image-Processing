%% Question 2
clc; clear all; close all;
tic;
%% Reading slices
slice_50 = im2double(imread("../data/slice_50.png"));
slice_50 = padarray(slice_50, [37, 19], 0, 'both');;
slice_51 = im2double(imread("../data/slice_51.png"));
slice_51 = padarray(slice_51, [37, 19], 0, 'both');
slice_52 = im2double(imread("../data/slice_52.png"));
slice_52 = padarray(slice_52, [37, 19], 0, 'both');
%% Creating measurements
random_angles = unifrnd(0, 180, 1, 36);
measurements_slice_50 = radon(slice_50, random_angles);
measurements_slice_51 = radon(slice_51, random_angles);
%% Part a - Filtered Back Projection using Ram-Lak Filter
[reconstructed_slice_50_ramlak, H] = iradon(measurements_slice_50, ...
    random_angles,'linear','Ram-Lak', 1, 255);
figure();
imshow(reconstructed_slice_50_ramlak);
colormap('gray');
title("Filtered Backprojection using Ram-Lak filter - Slice 50");
daspect([1 1 1]);
axis on;
axis tight;
colorbar;
% Slice 51
[reconstructed_slice_51_ramlak, H] = iradon(measurements_slice_51, ...
    random_angles,'linear','Ram-Lak', 1, 255);
figure();
imshow(reconstructed_slice_51_ramlak);
colormap('gray');
title("Filtered Backprojection using Ram-Lak filter - Slice 51");
daspect([1 1 1]);
axis on;
axis tight;
colorbar;
%% Part b - Independent CS-based reconstruction
addpath('./l1_ls_matlab');
y_50 = measurements_slice_50(:);
m = size(y_50, 1);
n = size(slice_50(:), 1);
projection_size = size(measurements_slice_50, 1);
A = forward_model(@idct2, projection_size, 255, random_angles);
At = forward_modelt(@dct2, projection_size, 255, random_angles);
lambda = 0.1;
rel_tol = 1e-6;
quiet = true;
[beta, status] = l1_ls(A, At, m, n, y_50, lambda,...
    rel_tol, quiet);
reconstructed_idpcs_slice_50 = idct2(reshape(beta, 255, 255));
figure();
imshow(reconstructed_idpcs_slice_50);
colormap('gray');
title("Independent CS-based Reconstruction - Slice 50");
daspect([1 1 1]);
axis on;
axis tight;
colorbar;
%
y_51 = measurements_slice_51(:);
m = size(y_51, 1);
n = size(slice_51(:), 1);
projection_size = size(measurements_slice_51, 1);
A = forward_model(@idct2, projection_size, 255, random_angles);
At = forward_modelt(@dct2, projection_size, 255, random_angles);
lambda = 0.1;
rel_tol = 1e-7;
quiet = true;
[beta, status] = l1_ls(A, At, m, n, y_50, lambda,...
    rel_tol, quiet);
reconstructed_idpcs_slice_51 = idct2(reshape(beta, 255, 255));
figure();
imshow(reconstructed_idpcs_slice_51);
colormap('gray');
title("Independent CS-based Reconstruction - Slice 51");
daspect([1 1 1]);
axis on;
axis tight;
colorbar;
%% Part c - Coupled Reconstruction
random_angles_set1 = unifrnd(0, 180, 1, 36);
random_angles_set2 = unifrnd(0, 180, 1, 36);
measurements_slice_50 = radon(slice_50, random_angles_set1);
measurements_slice_51 = radon(slice_51, random_angles_set2);
y_50 = measurements_slice_50(:);
y_51 = measurements_slice_51(:);
y_combined = [y_50; y_51];
m = size(y_combined, 1);
n = size(slice_50(:), 1) + size(slice_51(:), 1);
projection_size = size(measurements_slice_50, 1);
A = forward_model_coupled(@idct2, projection_size, 255, random_angles_set1, random_angles_set2);
At = forward_model_coupledt(@dct2, projection_size, 255, random_angles_set1, random_angles_set2);
lambda = 0.1;
rel_tol = 1e-6;
quiet = true;
[beta, status] = l1_ls(A, At, m, n, y_combined, lambda,...
    rel_tol, quiet);
beta1 = beta(1:0.5*n);
delta_beta1 = beta(1+0.5*n:end);
beta2 = beta1 + delta_beta1;
coupled_reconstructed_2_slice_50 = idct2(reshape(beta1, 255, 255));
coupled_reconstructed_2_slice_51 = idct2(reshape(beta2, 255, 255));
figure();
subplot(1, 2, 1),imshow(coupled_reconstructed_2_slice_50);
colormap('gray');
title("Coupled CS-based Reconstruction - Slice 50");
daspect([1 1 1]);
axis on;
axis tight;
colorbar;
subplot(1, 2, 2),imshow(coupled_reconstructed_2_slice_51);
colormap('gray');
title("Coupled CS-based Reconstruction - Slice 51");
daspect([1 1 1]);
axis on;
axis tight;
colorbar;
%% Reconstruction using 3 slices
random_angles_set1 = unifrnd(0, 180, 1, 36);
random_angles_set2 = unifrnd(0, 180, 1, 36);
random_angles_set3 = unifrnd(0, 180, 1, 36);
measurements_slice_50 = radon(slice_50, random_angles_set1);
measurements_slice_51 = radon(slice_51, random_angles_set2);
measurements_slice_52 = radon(slice_51, random_angles_set3);
y_50 = measurements_slice_50(:);
y_51 = measurements_slice_51(:);
y_52 = measurements_slice_52(:);
y_combined = [y_50; y_51; y_52];
m = size(y_combined, 1);
n = size(slice_50(:), 1) + size(slice_51(:), 1) + size(slice_52(:), 1);
projection_size = size(measurements_slice_50, 1);
A = forward_model_coupled_3(@idct2, projection_size, 255, random_angles_set1, random_angles_set2, random_angles_set3);
At = forward_model_coupled_3t(@dct2, projection_size, 255, random_angles_set1, random_angles_set2, random_angles_set3);
lambda = 0.1;
rel_tol = 1e-6;
quiet = true;
[beta, status] = l1_ls(A, At, m, n, y_combined, lambda,...
    rel_tol, quiet);
beta1 = beta(1:n/3);
delta_beta1 = beta(1+n/3:2*n/3);
delta_beta2 = beta(1+2*n/3:end);
beta2 = beta1 + delta_beta1;
beta3 = beta1 + delta_beta1 + delta_beta2;
coupled_reconstructed_3_slice_50 = idct2(reshape(beta1, 255, 255));
coupled_reconstructed_3_slice_51 = idct2(reshape(beta2, 255, 255));
coupled_reconstructed_3_slice_52 = idct2(reshape(beta3, 255, 255));
figure();
subplot(1, 3, 1),imshow(coupled_reconstructed_3_slice_50);
colormap('gray');
title("Coupled CS-based Reconstruction - Slice 50");
daspect([1 1 1]);
axis on;
axis tight;
colorbar;
subplot(1, 3, 2),imshow(coupled_reconstructed_3_slice_51);
colormap('gray');
title("Coupled CS-based Reconstruction - Slice 51");
daspect([1 1 1]);
axis on;
axis tight;
colorbar;
subplot(1, 3, 3),imshow(coupled_reconstructed_3_slice_52);
colormap('gray');
title("Coupled CS-based Reconstruction - Slice 52");
daspect([1 1 1]);
axis on;
axis tight;
colorbar;







%% toc;