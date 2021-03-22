classdef forward_model_coupled
    properties
        transform
        projection_size
        data_size
        projection_angles_1
        projection_angles_2
    end
    
    methods
        function obj = forward_model_coupled(transform,projection_size,...
                data_size, projection_angles_1, projection_angles_2)
            obj.transform = transform;
            obj.data_size = data_size;
            obj.projection_size = projection_size;
            obj.projection_angles_1 = projection_angles_1;
            obj.projection_angles_2 = projection_angles_2;
        end
        
        function product = mtimes(A, x)
            len_x = length(x);
            x1 = x(1:len_x/2);
            delta_x1 = x(0.5*len_x+1:end);
            x1 = reshape(x1, A.data_size, A.data_size);
            delta_x1 = reshape(delta_x1, A.data_size, A.data_size);
            beta1 = A.transform(x1);
            delta_beta1 = A.transform(delta_x1);
            R1U_beta1 = radon(beta1, A.projection_angles_1);
            R2U_beta1 = radon(beta1, A.projection_angles_2);
            R2U_delta_beta1 = radon(delta_beta1, A.projection_angles_2);
            product = [R1U_beta1(:); R2U_beta1(:) + R2U_delta_beta1(:)];
        end
    end
end