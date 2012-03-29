classdef MomentsToDiscrete
    % Converts a moments-based quantity to a discrete quantity.  However, since
    % we're limiting ourselves to isotropic scattering, we need only divide by
    % four pi.  We'll do this with a static vectorized function.
    
    properties (Access = private)
        %> Angular norm
        d_scale = 1.0 / 4.0 / pi; 
    end
    
    methods (Access = public)
        
        function this = MomentsToDiscrete(dim)
            if dim == 1
                this.d_scale = 0.5;
            end
        end

        function q = apply(this, v)
            q = v * this.d_scale;
        end
        
        function s = scale(this)
            s = 1/this.d_scale; 
        end
        
    end
    
    
end