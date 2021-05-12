% Returns the estimated gaussian model parameters
function [mu, sigma] = scs(y, seed, ratio, num_gaussians)
    num_iterations = 3;
    num_cols = size(y,2);
    [mu, sigma] = init(num_gaussians);
    model_set = containers.Map;

    % A map that contains set of signals which belong to ith gaussian model
    for i=1:num_gaussians
        model_set(int2str(i)) = zeros(1); % simple initialization, zeros will later be deleted
    end

    for iter=1:num_iterations
        % E-step
        rng(seed);
        for elem=1:num_cols
            phi = randn(floor(ratio*64), 64);
            prev_best_likelihood = -1e10;
            best_idx = 0;
            best_x_est = zeros(64, 1);
            for model_idx=1:num_gaussians
                sigma_j = reshape(sigma(:, model_idx), 64, 64);
                linear_decoder = sigma_j*phi'*inv(phi*sigma_j*phi');
                x_est(:, model_idx) = linear_decoder*y(:, elem)-linear_decoder*phi*mu(:, model_idx);
                log_likelihood = -log(det(sigma_j))-x_est(:, model_idx)'*inv(sigma_j)*x_est(:, model_idx);
                if prev_best_likelihood<log_likelihood
                    best_idx = model_idx;
                    best_x_est = x_est(:, model_idx)+mu(:, model_idx);
                    prev_best_likelihood = log_likelihood;
                end
            end
            model_set(int2str(best_idx)) = [model_set(int2str(best_idx)); best_x_est];
        end

        if iter==1  %remove the zero which was used for initialization
            for i=1:num_gaussians
                temp = model_set(int2str(i));
                temp = temp(2:end);
                model_set(int2str(i)) = temp;
            end
        end
        % M-step
        for model_idx=1:num_gaussians
            temp = model_set(int2str(model_idx));
            if size(temp, 2)>0
                temp = reshape(temp, 64, length(temp)/64);
                mu(:, model_idx) = mean(temp, 2);
                cov_matrix = cov(temp')+1e-3*eye(64);
                sigma(:, model_idx) = reshape(cov_matrix/(size(temp,2)), 64*64,1);
            end
        end
    end
end