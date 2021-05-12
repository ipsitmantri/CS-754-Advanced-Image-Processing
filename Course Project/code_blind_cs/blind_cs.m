clear all
clc
tic;

% change the name to lena, house or peppers to test other images
% for each iteration it takes almost 6mins and almost 40 mins
% for whole code to run
x = imread("../data/house.png"); 
x = imresize(x, 0.25); % image rescaled to 1/4 for faster reconstruction
x = double(rgb2gray(x));
ratio = [0.2, 0.25, 0.3, 0.4, 0.5, 0.6];

for r=1:length(ratio)
    seed = 1;
    [y, cnt] = sense(x, seed, ratio(r));   
    % y is the sensed signal and cnt is a 128*128 matrix where (i,j) element
    % is the number of times (i,j) pixel was reconstructed. Cnt will later
    % be used for averaging
    [theta dictionary] = reconstruct_bcs(y, seed, 15, 120);
    
    x_est = zeros(size(x));
    iter = 1;
    [rows, cols] = size(x);
    for i=1:rows-7
        for j=1:cols-7
            patch_est = dictionary*theta(:, iter);
            iter = iter + 1;
            x_est(i:i+7, j:j+7) = x_est(i:i+7, j:j+7) + reshape(patch_est, 8, 8);
        end
    end
    x_est = x_est./cnt;
    
    figure;
    imshow(mat2gray(x_est));
    title('Reconstructed Image');
    disp('rmse=');
    disp(norm(x_est(:)-x(:))/norm(x(:)));
    
    %change the saving names according to test image
    %as 'i' increases measurement ratio increases
    fullname = fullfile('../reconstructed_bcs_imgs', ['recons_house_' int2str(r) '.png']); 
    imwrite(mat2gray(x_est), fullname);
    
    % images saved as matfiles also for computing rmse
    save(['../mat_files/recons_house_bcs' int2str(r) '.mat'], 'x_est');
end

figure;
imshow(mat2gray(x));
title('Original Image')


toc;