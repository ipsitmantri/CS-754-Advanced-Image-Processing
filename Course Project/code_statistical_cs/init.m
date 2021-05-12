% EM-MAP initialization. Returns initial estimate of sigma and mu. 
% ith element in sigma is vector form of ith gaussian's sigma matrix estimate
% size(mu) = (64,19), size(sigma) = (64*64, 19)
function [mu, sigma] = init(num_gaussians)
    rng(5);
    sigma = zeros(64*64, num_gaussians);
    angles = pi*randn(num_gaussians, 1); %generates 19 random angles
    for l=1:num_gaussians
        theta = angles(l);
        X = ones(49, 49);
        for i=1:49
            for j=1:49
                if ((i-24)-(j-24)*tan(theta))<=1e-10
                    X(i,j)=0;
                end
            end
        end

        patch = zeros(64, 30);
        iter=1;
        for i=17:31
           patch_elem = X(i-3:i+4, 21:28);
           patch(:, iter) = patch_elem(:);
           iter = iter+1;
        end
        for j=17:31
           patch_elem = X(21:28, j-3:j+4);
           patch(:, iter) = patch_elem(:);
           iter = iter+1;
        end

        avg = mean(patch, 2);
        cov_matrix = zeros(64, 64);
        for i=1:30
            cov_matrix = cov_matrix + (patch(:, i)-avg)*(patch(:, i)-avg)';
        end
        cov_matrix = cov_matrix/(29)+1e-3*eye(64,64);
        sigma(:, l) = cov_matrix(:); % cov-matrix flattened and saved
    end
    mu = zeros(64, num_gaussians);
end