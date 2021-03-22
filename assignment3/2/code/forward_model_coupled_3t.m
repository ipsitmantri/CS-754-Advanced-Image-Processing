classdef forward_model_coupled_3t
    properties
        transform
        projection_size
        data_size
        projection_angles_1
        projection_angles_2
        projection_angles_3
    end
    
    methods
        function obj = forward_model_coupled_3t(transform, projection_size,...
                data_size, projection_angles_1, projection_angles_2, projection_angles_3)
            obj.transform = transform;
            obj.projection_size = projection_size;
            obj.data_size = data_size;
            obj.projection_angles_1 = projection_angles_1;
            obj.projection_angles_2 = projection_angles_2;
            obj.projection_angles_3 = projection_angles_3;
        end
        
        function product = mtimes(At, y)
            len_y = length(y);
            y1 = y(1:len_y/3);
            y2 = y(1+len_y/3:2*len_y/3);
            y3 = y(1+2*len_y/3:end);
            
            y1 = reshape(y1, At.projection_size, size(At.projection_angles_1, 2));
            y2 = reshape(y2, At.projection_size, size(At.projection_angles_2, 2));
            y3 = reshape(y3, At.projection_size, size(At.projection_angles_3, 2));
            
            x1 = iradon(y1, At.projection_angles_1, 'linear', ...
                'Ram-Lak', 1, At.data_size);
            d_x1 = iradon(y2, At.projection_angles_2, 'linear', ...
                'Ram-Lak', 1, At.data_size);
            d_x2 = iradon(y3, At.projection_angles_3, 'linear', ...
                'Ram-Lak', 1, At.data_size);
            
            beta1 = At.transform(x1);
            d_beta1 = At.transform(d_x1);
            d_beta2 = At.transform(d_x2);
            
            product = [beta1(:) + d_beta1(:) + d_beta2(:);...
                d_beta1(:);...
                d_beta2(:)];
        end
    end
end
        