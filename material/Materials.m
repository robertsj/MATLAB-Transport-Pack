classdef Materials < handle
    
    properties (Access = private)
        
        %> Number of groups
        d_number_groups = 0;
        
        %> Number of materials
        d_number_materials = 0;
        
        %> Downscatter switch (when true, upscatter ignored)
        d_downscatter = 0;
        
        % Total cross section [material, group]
        d_sigma_t = 0;
        
        %> Fission [material, group]
        d_nu_sigma_f = 0;
        
        %> Fission spectrum [material, group]
        d_chi = 0;
        
        %> Scatter [material, group<-, group')
        d_sigma_s = 0;
        
        %> Diffusion coefficient [material, group]
        d_diff_coef = 0;
        
        %> Scatter bounds applied to all materials [group, 2]
        d_scatter_bounds = 0;
        
        %> Upscatter cutoff.  Only groups equal to or above this cutoff are 
        %> subject to upscatter iterations.
        d_upscatter_cutoff = 1;
        
    end
    
    methods
       
        function obj = Materials(number_groups, ...
                                 number_materials, ...
                                 downscatter)
            if (number_groups < 1)
                error('Number of groups must be >= 1')
            end
            if (number_materials < 1)
                error('Number of materials must be >= 1')
            end
            if nargin < 3
                obj.d_downscatter = 0; % off by default
            else
                obj.d_downscatter = downscatter;
            end
            obj.d_number_groups = number_groups;
            obj.d_number_materials = number_materials;
            
            % Presize data arrays.
            obj.d_sigma_t    = zeros(number_materials, number_groups);
            obj.d_nu_sigma_f = zeros(number_materials, number_groups);
            obj.d_chi        = zeros(number_materials, number_groups);
            obj.d_sigma_s    = zeros(number_materials, number_groups, ...
                                     number_groups);
            obj.d_diff_coef  = zeros(number_materials, number_groups);                                 
                                 
            % Default: all groups are always added to the scatter source.
            obj.d_scatter_bounds = zeros(number_groups, 2);
            obj.d_scatter_bounds(:, 1) = 1;
            obj.d_scatter_bounds(:, 2) = number_groups;
            
        end
        
        % ----------------------------------------------------------------------
        % Public Interface
        % ----------------------------------------------------------------------
        
        % Setters
        
        function obj = set_sigma_t(obj, m, g, v)
            %check(obj, m, g)
            obj.d_sigma_t(m, g) = v;
        end
        
        function obj = set_nu_sigma_f(obj, m, g, v)
            %check(obj, m, g)
            obj.d_nu_sigma_f(m, g) = v;
        end
        
        function obj = set_chi(obj, m, g, v)
            %check(obj, m, g)
            obj.d_chi(m, g) = v;
        end
        
        function obj = set_sigma_s(obj, m, g, gp, v)
            %check(obj, m, g, gp)
            obj.d_sigma_s(m, g, gp) = v;
        end
        
        function obj = set_diff_coef(obj, m, g, v)
            %check(obj, m, g)
            obj.d_diff_coef(m, g) = v;
        end        
        
        % Setters (vectorized)
        
        function obj = set_sigma_t_v(obj, m, v)
            obj.d_sigma_t(m, :) = v;
        end
        
        function obj = set_nu_sigma_f_v(obj, m, v)
            obj.d_nu_sigma_f(m, :) = v;
        end
        
        function obj = set_chi_v(obj, m, v)
            obj.d_chi(m, :) = v;
        end
        
        function obj = set_sigma_s_v(obj, m, v)
            obj.d_sigma_s(m, :, :) = v;
        end   
        
        function obj = set_diff_coef_v(obj, m, v)
            obj.d_diff_coef(m, :) = v;
        end        
        
        % Getters
        
        function y = sigma_t(obj, m, g)
            %check(obj, m, g)
            y = obj.d_sigma_t(m, g);
        end
        
        function y = sigma_a(obj, m, g)
            %check(obj, m, g)
            y = obj.d_sigma_t(m, g) - sum(obj.d_sigma_s(m, :, g));
        end        
        
        function y = nu_sigma_f(obj, m, g)
            %check(obj, m, g)
            y = obj.d_nu_sigma_f(m, g);
        end
        
        function y = chi(obj, m, g)
            %check(obj, m, g)
            y = obj.d_chi(m, g);
        end
        
        function y = sigma_s(obj, m, g, gp)
            %check(obj, m, g, gp)
            y = obj.d_sigma_s(m, g, gp);
        end
        
        function y = diff_coef(obj, m, g)
            %check(obj, m, g)
            y = obj.d_diff_coef(m, g);
        end        
        
        function g = number_groups(obj)
           g = obj.d_number_groups; 
        end
        
        function b = lower(obj, g)
            b = obj.d_scatter_bounds(g, 1); 
        end
        
        function b = upper(obj, g)
            b = obj.d_scatter_bounds(g, 2); 
        end
        
        function b = downscatter(obj)
            b = obj.d_downscatter;
        end
        
        function b = upscatter_cutoff(obj)
            b = obj.d_upscatter_cutoff; 
        end
        
        function n = number_materials(obj)
           n = obj.d_number_materials; 
        end
        
        function obj = finalize(obj)
            
            % Set the scatter group bounds.  For each group, we compute the
            % lowest index (highest energy) that leads to downscatter.  We also
            % compute the highest index (lowest energy) that upscatters into the
            % group.  Knowing these bounds eliminates a bit of computation in
            % computing the scattering source.
            for g = 1:obj.d_number_groups
                lower = g;
                upper = g;
                
                for m = 1:obj.d_number_materials
                    
                    % Downscatter from gp to g
                    for gp = 1:g
                        if obj.d_sigma_s(m, g, gp) > 0
                            lower = min(gp, lower);
                        end
                    end
                    
                    % Upscatter from gp to g
                    for gp = g:obj.d_number_groups
                        if obj.d_sigma_s(m, g, gp) > 0
                            upper = max(gp, upper);
                        end
                    end             
                    
                end
                obj.d_scatter_bounds(g, 1) = lower; 
                obj.d_scatter_bounds(g, 2) = upper; 
            end
            
            % Go through the scatter bounds for each g.  If for some g, the
            % upper scatter bound is larger than g, then upscatter exists
            % from that lower energy group.  The first group g for which
            % this occurs is the upscatter cutoff.
            obj.d_upscatter_cutoff = obj.d_number_groups + 1;
            for g = 1:obj.d_number_groups
                if obj.d_scatter_bounds(g, 2) > g
                    obj.d_upscatter_cutoff = g;
                    break;
                end
            end
            
            % If our materials have no upscatter, then we set the
            % downscatter-only flag.
            if obj.d_upscatter_cutoff == obj.d_number_groups + 1
                if obj.d_downscatter == 0
                   warning('user:input', ...
                       ['Upscatter is being turned off since', ...
                        ' no upscatter exists in the data.']);
                end
                obj.d_downscatter = 1;
            end
            
        end
        
        %> @brief  Print all material data in a nice table.
        function print_materials(obj)
            
            
        end
  
        
    end
    
    methods (Access = private)
        
        function check(obj, m, g, gp)
            if (m < 1 || m > obj.d_number_materials)
                error('Invalid material index');
            end
            if (g < 1 || g > obj.d_number_groups)
                error('Invalid group index');
            end
            if (nargin == 4)
                if (gp < 1 || gp > obj.d_number_groups)
                    error('Invalid group prime index');
                end
            end
        end
        
    end % end private methods
    
    
end