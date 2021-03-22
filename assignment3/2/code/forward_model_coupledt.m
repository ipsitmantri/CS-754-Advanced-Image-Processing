classdef forward_model_coupledt
    properties
        transform
        projection_size
        data_size
        projection_angles_1
        projection_angles_2
    end
    
    methods
        function obj = forward_model_coupledt(transform, projection_size,...
                data_size, projection_angles_1, projection_angles_2)
            obj.transform = transform;
            obj.projection_size = projection_size;
            obj.data_size = data_size;
            obj.projection_angles_1 = projection_angles_1;
            obj.projection_angles_2 = projection_angles_2;
        end
        
        function product = mtimes(At, y)
            len_y = length(y);
            y1 = y(1:len_y/2);
            y2 = y(0.5*len_y+1:end);
            
            y1 = reshape(y1, At.projection_size, size(At.projection_angles_1,2));
            y2 = reshape(y2, At.projection_size, size(At.projection_angles_2,2));
            
            beta1 = iradon(y1, At.projection_angles_1, 'linear',...
                'Ram-Lak', 1, At.data_size);
            delta_beta1 = iradon(y2, At.projection_angles_2, 'linear',...
                'Ram-Lak', 1, At.data_size);
            
            x1 = At.transform(beta1);
            delta_x1 = At.transform(delta_beta1);
            
            product = [x1(:) + delta_x1(:); delta_x1(:)];
        end
    end
end
            
        