x = imread("barbara256.png");
x = double(x);
[rows cols] = size(x);
padded_x = padarray(x,[7, 7],0,"both"); %padding x for easier avergaing
myNumOfColors = 200;
myColorScale = [ [0:1/(myNumOfColors-1):1]' ,[0:1/(myNumOfColors-1):1]' , [0:1/(myNumOfColors-1):1]' ];

dct_2d = kron(dctmtx(8)',dctmtx(8)');
stride = 1;
[rows cols] = size(padded_x);
rng(4);
phi = randn(32, 64);
x_est = zeros(size(padded_x));
A = phi*dct_2d;
max_eigen = max(eig(A'*A));
alpha = max_eigen+1;
for i=1:stride:rows-7
    for j=1:stride:cols-7
        patch = padded_x(i:i+7, j:j+7);
        y = phi*patch(:);
        lambda = 1;
        theta_est = zeros(64,1);
        % Ista algo
        for iter=1:200
            theta_est = soft(theta_est+(1/alpha)*A'*(y-A*theta_est), lambda/(2*alpha)); %adding the newly estimated pixels with previous one
        end
        est_patch = dct_2d*theta_est;
        est_patch = reshape(est_patch, 8, 8);
        x_est(i:i+7, j:j+7) = x_est(i:i+7, j:j+7) + est_patch;
    end
end   
x_est = x_est(8:263, 8:263)/64; %because of padding averaging is same for all
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