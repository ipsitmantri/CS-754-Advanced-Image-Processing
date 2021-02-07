function I_reconstructed = omp_video_reconstruction(I, S, patch_size, T, epsilon)
    p = patch_size;
    dct_basis = dctmtx(p);
    dct_basis_2D = kron(dct_basis, dct_basis);
    psi = kron(eye(T), dct_basis_2D);
    I_reconstructed = zeros(size(I,1), size(I,2), T);
    counts = zeros(size(I,1), size(I,2), T);
    for i=1:size(I,1) - p + 1
        for j=1:size(I,2) - p + 1
            patch = I(i:i+p-1, j:j+p-1);
            patch_S = S(i:i+p-1, j:j+p-1,:);
            phi = [];
            for t=1:size(S, 3)
                patch_St = patch_S(:, :, t);
                phi_t = diag(reshape(patch_St, p*p, 1));
                phi = [phi phi_t];
            end
            A = phi * psi';
            y = reshape(patch, p*p, 1);
            [theta, ~] = omp(y, A, epsilon);
            i_t = psi' * theta;
            for k=1:T
                f = reshape(i_t((k-1)*p*p+1:k*p*p,:), p, p);
                I_reconstructed(i:i+p-1, j:j+p-1,k) = I_reconstructed(i:i+p-1, j:j+p-1,k) + f;
                counts(i:i+p-1, j:j+p-1,k) = counts(i:i+p-1, j:j+p-1,k) + ones(p,p);
            end
%            disp([i, j]);
        end
    end
    I_reconstructed = I_reconstructed ./ counts;
end
function [theta_opt, T] = omp(y, A, epsilon)
    % Orthogonal Matching Pursuit
    r = y;
    theta = zeros(size(A, 2), 1);
    i = 0;
    T = [];
    A_N = normc(A);
    while (norm(r) > epsilon & i < 64)
        [~, j] = max(abs(r' * A_N));
        T = [T j];
        i = i + 1;
        A_T_i = A(:, T);
        theta(T) = pinv(A_T_i) * y;
        r = y - A_T_i * theta(T);
    end
    theta_opt = theta;
end