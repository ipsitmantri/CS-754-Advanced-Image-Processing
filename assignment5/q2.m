clear all
clc
n = 128;
rng(0);
U = RandOrthMat(n);   %function taken from internet
m = [40, 50, 64, 80, 100, 120];
alpha = [0, 1, 2, 3, 4];
figure;
for l=1:length(alpha)
    diagonal = zeros(128,1);
    for i=1:128
        diagonal(i) = sqrt(i^(-alpha(l)));
    end
    A = U*diag(diagonal);
    rng(1);
    rmse = zeros(length(m),length(alpha));
    for i = 1:length(m)
        for j=1:10
            x = A*randn(128,1);    %Generates random numbers from the distribution N(0, AA')
            phi = sqrt(1/m(i))*randn(m(i), 128);   %sensing matrix
            sensed_x = phi*x;
            abs_avg = mean(abs(sensed_x));
            sigma = 0.01*abs_avg;
            noise = sigma*randn(m(i),1);
            y = sensed_x + noise;
            recons_x = (inv(phi'*phi+sigma^2*(U*diag(1./(diagonal.^2))*U')))*phi'*y;
            rmse(i, l) = rmse(i, l)+norm(recons_x-x)/sqrt(128);
        end
        rmse(i, l) = rmse(i, l)/10;
    end
    plot(m, log(rmse(:, l)+0.01), '-o');   %log scale to better visulaize rmse variation
    hold on
end
ylabel('log(RMSE+0.01)');
xlabel('# of measurements(m)')
legend('alpha=0', 'alpha=1', 'alpha=2', 'alpha=3', 'alpha=4');
