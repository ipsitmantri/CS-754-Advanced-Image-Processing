function x_est = reconstruct_CS(y, cnt, psi, seed, ratio, rows, cols)
    rng(seed);
    it = 1;
    x_est = zeros(rows, cols);
    x_vec = zeros(64, size(y, 2));
    for i=1:size(y, 2)
%         disp(i);
        phi = randn(floor(ratio*64), 64);
        A = phi * psi;
        max_eigen = max(eig(A'*A));
        alpha = max_eigen+1;
        theta_est = zeros(225,1);
        lambda = 100;
        % Ista algo
        for iter=1:100
            theta_est = soft(theta_est+(1/alpha)*A'*(y(:,i)-A*theta_est), lambda/(2*alpha));
        end
%         [thetai, status] = l1_ls(A, y(:,i), 0.01, 0.01, true);
        x_vec(:, i) = psi * theta_est;
    end
    
    for p=1:rows-7
        for q=1:cols-7
            x_est(p:p+7, q:q+7) = x_est(p:p+7, q:q+7) + reshape(x_vec(:, it), 8, 8);
            it = it + 1;
        end
    end
    x_est = max(x_est./cnt, 0);
end
