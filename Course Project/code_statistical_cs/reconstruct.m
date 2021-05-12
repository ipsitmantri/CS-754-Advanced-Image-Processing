% returns the estimated signal using the given gaussian model parameters
% reconstructs using linear decoder and MAP estimation
function x_est = reconstruct(y, seed, ratio, num_gaussians, mu, sigma, rows, cols)
    iter = 1;
    rng(seed);
    x_est = zeros(rows, cols);
    for i=1:rows-7
        for j=1:cols-7
            phi = randn(floor(ratio*64), 64);
            prev_best_likelihood = -1e10;
            best_idx = 0;
            best_x_est = zeros(64, 1);
            for model_idx=1:num_gaussians
                sigma_j = reshape(sigma(:, model_idx), 64, 64);
                linear_decoder = sigma_j*phi'*inv(phi*sigma_j*phi');
                patch = linear_decoder*y(:, iter)-linear_decoder*phi*mu(:, model_idx);
                log_likelihood = -log(det(sigma_j))-patch'*inv(sigma_j)*patch;
                if prev_best_likelihood<log_likelihood
                    best_idx = model_idx;
                    best_x_est = patch+mu(:, model_idx);
                    prev_best_likelihood = log_likelihood;
                end
            end
            x_est(i:i+7, j:j+7) = x_est(i:i+7, j:j+7) + reshape(best_x_est, 8, 8);
            iter=iter+1;
        end
    end
end