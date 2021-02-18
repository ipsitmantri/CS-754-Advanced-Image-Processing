% ista for haar wavelet basis
% We need to evaluate AT*y-AT*A*theta, where A = phi*idwt2
function theta_est = ista_dwt(y, phi, alpha, lambda)
    theta_est = zeros(64,1);
    for iter=1:20
        [A, B, C, D] = dwt2(reshape(phi'*y, 8, 8), 'db1');
        ATy = [A(:); B(:); C(:); D(:)]; 
        A_est = reshape(theta_est(1:16), 4, 4);
        B_est = reshape(theta_est(17:32), 4, 4);
        C_est = reshape(theta_est(33:48), 4, 4);
        D_est = reshape(theta_est(49:64), 4, 4);
        im = idwt2(A_est,B_est,C_est,D_est,'db1');
        M = phi'*phi*im(:);
        [A, B, C, D] = dwt2(reshape(M, 8, 8), 'db1');
        ATAtheta = [A(:); B(:); C(:); D(:)];
        theta_est = soft(theta_est+(1/alpha)*(ATy - ATAtheta), lambda/(2*alpha));
    end
end