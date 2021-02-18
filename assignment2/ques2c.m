%% This code takes more than an hour to run.
x = imread("barbara256.png");
x = double(x);
[rows cols] = size(x);
padded_x = padarray(x,[7, 7],0,"both"); %padding x for easier avergaing
myNumOfColors = 200;
myColorScale = [ [0:1/(myNumOfColors-1):1]' ,[0:1/(myNumOfColors-1):1]' , [0:1/(myNumOfColors-1):1]' ];

stride = 1;
[rows cols] = size(padded_x);
rng(4);
phi = randn(32, 64);
x_est = zeros(size(padded_x));
alpha = 500;    %choosen by trail 
for i=1:stride:rows-7
    for j=1:stride:cols-7
        patch = padded_x(i:i+7, j:j+7);
        y = phi*patch(:);
        lambda = 1;
        theta_est = ista_dwt(y, phi, alpha, lambda); %needed to reshape a lot of matrix hence improvised the ista algorithm
        A_est = reshape(theta_est(1:16), 4, 4);
        B_est = reshape(theta_est(17:32), 4, 4);
        C_est = reshape(theta_est(33:48), 4, 4);
        D_est = reshape(theta_est(49:64), 4, 4);
        est_patch = idwt2(A_est,B_est,C_est,D_est,'db1'); %db1 == Haar wavelet basis
        x_est(i:i+7, j:j+7) = x_est(i:i+7, j:j+7) + est_patch;
    end
end   
x_est = x_est(8:263, 8:263)/64;
rmse = norm(x-x_est, 'fro')/norm(x, 'fro')
%%
figure(1), imagesc (single (x)); 
title('Original Image')
colormap (myColorScale);
%colormap jet;
daspect ([1 1 1]);
axis tight;
colorbar

figure(2), imagesc (single (x_est)); 
title('Reconstructed Image')
colormap (myColorScale);
%colormap jet; 
daspect ([1 1 1]);
axis tight;
colorbar