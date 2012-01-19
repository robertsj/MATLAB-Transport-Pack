classdef MomentsToDiscrete
    % Converts a moments-based quantity to a discrete quantity.  However, since
    % we're limiting ourselves to isotropic scattering, we need only divide by
    % four pi.  We'll do this with a static vectorized function.
    
    properties (Access = private)
        %> Angular norm
        d_scale = 1.0 / 4.0 / pi; 
    end
    
    methods (Access = public)
        
        function obj = MomentsToDiscrete(dim)
            if dim == 1
                obj.d_scale = 0.5;
            end
        end

        function q = apply(obj, v)
            q = v * obj.d_scale;
        end
        
    end
    
    
end