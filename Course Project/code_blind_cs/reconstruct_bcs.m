% this function reconstructs image as well as the dictionary from the
% sensed signal using the blind compressed sensing algorithm
function  [theta, dic] = reconstruct_bcs(y, seed, s, atoms)
    [rows, cols] = size(y);
    dic = randn(64, atoms);
    theta = zeros(atoms, cols);
    for iter = 1:5
        % sparse code
        rng(1);
        for i = 1:cols
            phi = randn(rows, 64);
            theta(:, i) = omp(phi*dic, y(:, i), s);
        end
        % dictionary atoms update
        for k = 1:atoms
            rng(1);
            term1 = zeros(64, 64);
            term2 = zeros(64, 1);
            y_res = zeros(size(y));
            for i = 1:cols
                phi = randn(rows, 64);
                y_res(:, i) = y(:, i) - phi*dic*theta(:, i) + phi*dic(:, k)*theta(k,i);
                term1 = term1 + phi'*phi*theta(k, i)^2;
                term2 = term2 + phi'*y_res(:, i)*theta(k, i);
            end
            dic(:, k) = inv(term1)*term2;
            dic(:, k) = dic(:, k)/norm(dic(:, k));
         end
    end
    rng(1);
    for i = 1:cols
        phi = randn(rows, 64);
        theta(:, i) = omp(phi*dic, y(:, i), s);
    end
end