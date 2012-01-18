classdef MomentsToDiscrete
    % Converts a moments-based quantity to a discrete quantity.  However, since
    % we're limiting ourselves to isotropic scattering, we need only divide by
    % four pi.  We'll do this with a static vectorized function.
    
    properties (Constant)
        
       d_scale = 1.0 / 4.0 / pi; 
       
    end
    
    methods (Static)
        
        function q = apply(v)
            q = v * MomentsToDiscrete.d_scale;
        end
        
    end
    
    
end