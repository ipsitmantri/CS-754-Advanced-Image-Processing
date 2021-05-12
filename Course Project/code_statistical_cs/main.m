clear all
clc
tic;

% change the name to lena, house or peppers to test other images
% for each iteration it takes almost 6mins and almost 40 mins
% for whole code to run
x = imread("../data/lake.png"); 
x = imresize(x, 0.25); % image rescaled to 1/4 for faster reconstruction
x = double(rgb2gray(x));
ratio = [0.2, 0.25, 0.3, 0.4, 0.5, 0.6];

for i=1:length(ratio)
    seed = 1;
    [y, cnt] = sense(x, seed, ratio(i));   
    %y is the sensed signal and cnt is a 128*128 matrix where (i,j) element
    %is the number of times (i,j) pixel was reconstructed. Cnt will later
    %be used for averaging

    num_gaussians = 19;
    [mu, sigma] = scs(y, seed, ratio(i), num_gaussians); %optimal gaussian parameters returned

    x_est = reconstruct(y, seed, ratio(i), num_gaussians, mu, sigma, size(x,1), size(x,2));
    x_est = max(x_est./cnt, 0);  %pixels which are less than 0 are clipped to zero
    
    figure;
    imshow(mat2gray(x_est));
    title('Reconstructed Image');
    
    %change the saving names according to test image
    %as 'i' increases measurement ratio increases
    fullname = fullfile('../reconstructed_scs_imgs', ['recons_lake_' int2str(i) '.png']); 
    imwrite(mat2gray(x_est), fullname);
    
    % images saved as matfiles also for computing rmse
    save(['../mat_files/recons_lake_' int2str(i) '.mat'], 'x_est');
end

figure;
imshow(mat2gray(x));
title('Original Image')


toc;