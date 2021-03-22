classdef forward_model
    properties
        transform
        data_size
        projection_size
        projection_angles
    end
    
    methods
        function obj = forward_model(transform, measurement_size,...
                data_size, projection_angles)
            obj.transform = transform;
            obj.projection_size = measurement_size;
            obj.data_size = data_size;
            obj.projection_angles = projection_angles;
        end
        
        function product = mtimes(A, x)
            x = reshape(x, A.data_size, A.data_size);
            beta = A.transform(x);
            product = radon(beta, A.projection_angles);
            product = product(:);
        end
    end
end
            
            