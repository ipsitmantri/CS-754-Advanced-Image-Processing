%% Clearing Workspace and Loading Toolboxes
clc; clear all; close all;
addpath('../code_statistical_cs/');
addpath('../ompbox/');
addpath('../ksvdbox/');
tic;
%% Combining the 200 train and 100 test images from Berkely Segmentation Dataset
train_images = dir(fullfile('../BSDS300/images/train/', '*.jpg'));
test_images = dir(fullfile('../BSDS300/images/test/', '*.jpg'));
image_names = {};
for i=1:size(train_images, 1)
    image_names = [image_names, "../BSDS300/images/train/" + convertCharsToStrings(train_images(i).name)];
end
for i=1:size(test_images, 1)
    image_names = [image_names, "../BSDS300/images/test/" + convertCharsToStrings(test_images(i).name)];
end
%% Preparing the training data to apply K-SVD
Y = []; % a matrix of size 64 x 720,000
for i=1:size(image_names, 2)
    img = rgb2gray(imread(image_names(i)));
    fun = @(block_struct) block_struct.data(:); 
    Yi = blockproc(img(2:end, 2:end), [8, 8], fun, 'PadPartialBlocks', true);
    Yi = reshape(Yi, 64, []);
    Y = [Y, Yi];
end
%%
n = size(Y, 1); % no. of rows in dictionary
m = 225; % no. of columns in dictionary
L = size(Y, 2); % no. of training samples
k = 3; % sparsity level of learned signal representation
params.data = double(Y); % training data
params.Tdata = k;
params.dictsize = m;
params.iternum = 50; % no. of iterations for K-SVD
params.memusage = 'high';
[Dksvd,g,err] = ksvd(params,''); % Using the KSVD Toolbox

figure; plot(err); title('K-SVD error convergence');
xlabel('Iteration'); ylabel('RMSE');

fprintf("  Dictionary size: %d x %d", n, m);
fprintf("  Number of examples: %d", L);
%%
seed = 1;
rng(seed);
psi = Dksvd;
ratio = [0.2, 0.25, 0.3, 0.4, 0.5, 0.6]; % fration of measurements taken from 8x8 patches
x = imread('../data/peppers.png'); % test image
x = imresize(x, 0.25); % resizing it to increase computational speed
x = double(rgb2gray(x));
rmse = []; % stores RMSE for different sensing matrices
for i=1:size(ratio, 2)
    [y, cnt] = sense(x, seed, ratio(i)); % sensing measurements
    x_est = reconstruct_CS(y, cnt, psi, seed, ratio(i), size(x, 1), size(x, 2)); % CS reconstruction using the learned overcomplete dictionary
                                                                                 % and ISTA       
    
    figure();
    imshow(mat2gray(x_est));
    title('Reconstructed Image');
    
    %change the saving names according to test image
    %as 'i' increases measurement ratio increases
    fullname = fullfile('../reconstructed_cs_imgs', ['recons_peppers_' int2str(i) '.png']); 
    imwrite(mat2gray(x_est), fullname);
    
    % images saved as matfiles also for computing rmse
    save(['../mat_files/recons_peppers_cs_' int2str(i) '.mat'], 'x_est');
    rmse = [rmse, norm(x_est(:)-x(:))/norm(x(:))];
end
figure;
imshow(mat2gray(x));
title('Original Image');
figure();
plot(ratio, rmse, 'bo-');
xlabel('ratio');
ylabel('RMSE');
title('RMSE vs ratio for Peppers');
%%
toc;