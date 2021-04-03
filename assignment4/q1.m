clc
clear all
%% part a and b
s = [5:5:20, 30:10:100];   %sparsity level
sigma_amp = [0:0.1:0.4, 0.6:0.3:2];  %noise standard deviation

basis1 = dctmtx(256)';      %dct basis 
basis2 = eye(256);      %temporal basis  
noise = zeros(256,1);
A = [basis1, basis2];
max_eigen = max(eig(A'*A));
alpha = max_eigen+0.1;
rmse1 = zeros(length(s),length(sigma_amp));
rmse2 = zeros(length(s),length(sigma_amp));

for i=1:length(s)
    rng(0);
    p = randperm(256,s(i));        %select any s(i) frequencies from 256 indices
    rng(0);
    freq = zeros(256,1);        %frequency vector
    amplitude = randi([5 15],1,s(i));
    freq(p) = amplitude;
    f1 = basis1*freq;

    f2 = zeros(256,1);
    rng(1);
    p = randperm(256,s(i));   %select any s(i) time indices from 256 indices
    rng(1);
    values = randi([5 15],1,s(i));
    f2(p) = values;
    
    for j=1:length(sigma_amp)
        rng(2);
        noise = randn(256,1)*sigma_amp(j)*(mean(f1)+mean(f2));
        f = f1+f2+noise;
        J = zeros(200,1); % cost function 
        x_est = zeros(512,1);
        lambda = 1;   %regularization parameter
        for iter=1:200
            x_est = soft(x_est+A'*(f-A*x_est)./alpha,lambda/(2*alpha));  %ISTA
            J(iter) = sum((f-A*x_est).^2) + lambda*sum(abs(x_est));  %cost function
        end
        %figure, plot(1:200, J);
        rmse1(i, j) = norm(f1-basis1*x_est(1:256))/sqrt(length(f));
        rmse2(i, j) = norm(f2-x_est(257:end))/sqrt(length(f));
    end
end
figure;
for i=1:length(sigma_amp)
    plot(s, rmse1(:, i), 'o-');
    hold on;
end
xlabel('sparsity(s)');
ylabel('RMSE');
Legend=cell(length(sigma_amp),1);
for i=1:length(sigma_amp)
    Legend{i}=strcat('std factor=', num2str(sigma_amp(i)));
end
legend(Legend)

figure;
for i=1:length(sigma_amp)
    plot(s, rmse2(:, i), 'o-');
    hold on;
end
xlabel('sparsity(s)');
ylabel('RMSE');
Legend=cell(length(sigma_amp),1);
for i=1:length(sigma_amp)
    Legend{i}=strcat('std factor=', num2str(sigma_amp(i)));
end
legend(Legend)
%% partc s=20, sigma_amp=0.5
s = 20;   %sparsity level
sigma_amp = 0.1;  %noise standard deviation
k = [0.1:0.1:1, 2:10];

basis1 = dctmtx(256)';      %dct basis 
basis2 = eye(256);      %temporal basis  
noise = zeros(256,1);
A = [basis1, basis2];
max_eigen = max(eig(A'*A));
alpha = max_eigen+0.1;
rmse1 = zeros(length(k),1);
rmse2 = zeros(length(k),1);

for i=1:length(k)
    rng(0);
    p = randperm(256,s);        %select any s frequencies from 256 indices
    rng(0);
    freq = zeros(256,1);        %frequency vector
    amplitude = randi([5 15],1,s);
    freq(p) = amplitude;
    f1 = basis1*freq;

    f2 = zeros(256,1);
    rng(1);
    p = randperm(256,s);   %select any s time indices from 256 indices
    rng(1);
    values = k(i)*amplitude; %magnitude of f2 is k*f1
    f2(p) = values;
    
    rng(2);
    noise = randn(256,1)*sigma_amp*(mean(f1)+mean(f2));
    f = f1+f2+noise;
    J = zeros(200,1); % cost function 
    x_est = zeros(512,1);
    lambda = 1.75;   %regularization parameter
    for iter=1:200
        x_est = soft(x_est+A'*(f-A*x_est)./alpha,lambda/(2*alpha));  %ISTA
        J(iter) = sum((f-A*x_est).^2) + lambda*sum(abs(x_est));  %cost function
    end
    %figure, plot(1:200, J);
    rmse1(i) = norm(f1-basis1*x_est(1:256))/sqrt(length(f));
    rmse2(i) = norm(f2-x_est(257:end))/sqrt(length(f));
end
figure;
plot(k, rmse1, 'o-');
hold on;
plot(k, rmse2, 'o-');
xlabel('Magtindue ratio');
ylabel('RMSE');
legend('f1', 'f2');

