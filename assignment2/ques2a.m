x = imread("barbara256.png");
x = double(x);
padded_x = padarray(x,[7, 7],0,"both"); %padding x for easier avergaing
myNumOfColors = 200;
myColorScale = [ [0:1/(myNumOfColors-1):1]' ,[0:1/(myNumOfColors-1):1]' , [0:1/(myNumOfColors-1):1]' ];

rng(0);
n = 2*randn(size(padded_x));
y = padded_x+n; %noise with variance 4

dct_2d = kron(dctmtx(8)',dctmtx(8)');
stride = 1; 
[rows cols] = size(padded_x);
phi = eye(64);  %measurement matrix
x_est = zeros(size(padded_x));
A = phi*dct_2d;
max_eigen = max(eig(A'*A));
alpha = max_eigen+1;
for i=1:stride:rows-7
    for j=1:stride:cols-7
        y_patch = y(i:i+7, j:j+7);
        y_patch = y_patch(:);
        lambda = 1;
        theta_est = zeros(64,1);
        % Ista algo
        for iter=1:100
            theta_est = soft(theta_est+(1/alpha)*A'*(y_patch-A*theta_est), lambda/(2*alpha));
        end
        est_patch = dct_2d*theta_est;
        est_patch = reshape(est_patch, 8, 8);
        x_est(i:i+7, j:j+7) = x_est(i:i+7, j:j+7) + est_patch; %adding the newly estimated pixels with previous one
    end
end   

x_est = x_est(8:263, 8:263)/64;  %because of padding averaging is same for all
rmse = norm(x-x_est, 'fro')/norm(x, 'fro')
error = norm(x-y(8:263, 8:263), 'fro')/norm(x, 'fro')

%%
figure;
subplot(1,3,1), imagesc (single (x)); 
title('Original Image')
colormap (myColorScale);
%colormap jet;
daspect ([1 1 1]);
axis tight;
colorbar

subplot(1,3,2), imagesc (single (y)); 
title('Noisy Image')
colormap (myColorScale);
%colormap jet;
daspect ([1 1 1]);
axis tight;
colorbar

subplot(1,3,3), imagesc (single (x_est)); 
title('Reconstructed Image')
colormap (myColorScale);
%colormap jet; 
daspect ([1 1 1]);
axis tight;
colorbar