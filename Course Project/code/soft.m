%soft function
function theta_new = soft(y, lambda)
    theta_new = zeros(size(y));
    for i=1:length(y)
        if y(i)>=lambda
            theta_new(i) = y(i)-lambda;
        elseif y(i)<=-lambda
            theta_new(i) = y(i)+lambda;
        else
            theta_new(i) = 0;
        end
end
