classdef forward_modelt
    properties
        transform
        data_size
        projection_size
        projection_angles
    end
    
    methods
        function obj = forward_modelt(transform, measurement_size,...
                data_size, projection_angles)
            obj.transform = transform;
            obj.projection_size = measurement_size;
            obj.data_size = data_size;
            obj.projection_angles = projection_angles;
        end
        
        function product = mtimes(At, y)
            y = reshape(y, At.projection_size, size(At.projection_angles,2));
            beta = iradon(y, At.projection_angles, 'linear',...
'Ram-Lak', 1, At.data_size);
            x = At.transform(beta);
            product = x(:);
        end
    end
end
            
            