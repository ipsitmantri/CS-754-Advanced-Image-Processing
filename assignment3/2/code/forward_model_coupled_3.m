classdef forward_model_coupled_3
    properties
        transform
        projection_size
        data_size
        projection_angles_1
        projection_angles_2
        projection_angles_3
    end
    
    methods
        function obj = forward_model_coupled_3(transform, projection_size, data_size, projection_angles_1, projection_angles_2, projection_angles_3)
            obj.transform = transform;
            obj.projection_size = projection_size;
            obj.data_size = data_size;
            obj.projection_angles_1 = projection_angles_1;
            obj.projection_angles_2 = projection_angles_2;
            obj.projection_angles_3 = projection_angles_3;
        end
        
        function product = mtimes(A, beta)
            len_b = length(beta);
            beta1 = beta(1:len_b/3);
            d_beta1 = beta(len_b/3+1:2*len_b/3);
            d_beta2 = beta(2*len_b/3+1:end);
            beta1 = reshape(beta1, A.data_size, A.data_size);
            d_beta1 = reshape(d_beta1, A.data_size, A.data_size);
            d_beta2 = reshape(d_beta2, A.data_size, A.data_size);
            
            x1 = A.transform(beta1);
            d_x1 = A.transform(d_beta1);
            d_x2 = A.transform(d_beta2);
            
            R1x1 = radon(x1, A.projection_angles_1);
            R2x1 = radon(x1, A.projection_angles_2);
            R3x1 = radon(x1, A.projection_angles_3);
            R2d_x1 = radon(d_x1, A.projection_angles_2);
            R3d_x2 = radon(d_x2, A.projection_angles_3);
            
            product = [R1x1(:); R2x1(:) + R2d_x1(:); R3x1(:) ...
                + R3d_x2(:)];
        end
    end
end
            
        