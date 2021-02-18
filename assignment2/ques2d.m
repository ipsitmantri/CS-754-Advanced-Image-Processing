%given constants
rng(6);
indx = randperm(100);
h = [1; 2; 3; 4; 3; 2; 1]/16;

%required constants
indx = indx(51:60);
x = zeros(100, 1);
rng(2);
x(indx) = 5*rand(10,1)+1;
energy = sum(x(indx).^2);
sigma = 0.05*sqrt(energy);
rng(1);
n = randn(106, 1)*sigma;
y = conv(x, h)+n;
A = zeros(106, 106);
cols = 7;
for row=1:106
    A(row, cols-6:cols) = h;
    cols=cols+1;
end
A = A(:, 7:106);

%Ista Algorithm
lambda = 0.75;
x_est = zeros(100,1);
max_eigen = max(eig(A'*A));
alpha = max_eigen+0.1;
J = zeros(400,1); % cost function 
for iter=1:400
    x_est = soft(x_est+A'*(y-A*x_est)./alpha,lambda/(2*alpha));
    J(iter) = sum((y-A*x_est).^2) + lambda*sum(abs(x_est));
end
rmse = norm(x-x_est)/norm(x)
figure, plot(1:106, y);
figure, plot(1:400, J);
figure, plot(1:100, x_est);
figure, plot(1:100, x);