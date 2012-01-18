%> @file  Fixed.m
%> @brief Fixed class definition.
% ==============================================================================
%> @brief  Fixed source multigroup solver.
%
%> The multigroup transport equations are solved via Gauss-Seidel with
%> iterations over the upscatter block. 
% ==============================================================================
classdef Fixed
    
    properties
        d_input                 % User input.
        d_mesh                  % Cartesian mesh.
        d_mat                   % Materials.
        d_quadrature            % Angular mesh.
        d_source                % External volume source.
        d_boundary_source       % External boundary source.
        d_fission_source        % Fission source.
        d_state                 % State variables.
        d_inner_solver 
        d_boundary
    end
    
    methods
        
        % ======================================================================
        %> @brief Class constructor
        %
        %> @param input             Number of fine mesh per x coarse mesh.
        %> @param state             Number of fine mesh per x coarse mesh.
        %> @param boundary          Number of fine mesh per x coarse mesh.
        %> @param mesh              Number of fine mesh per y coarse mesh.
        %> @param mat               Coarse mesh edges along x axis.
        %> @param quadrature        Coarse mesh edges along y axis.
        %> @param external_source 	Coarse mesh material map.
        %> @param fission_source 	Coarse mesh material map.
        %>
        %> @return Instance of the Fixed class.
        % ======================================================================
        function obj = Fixed(input,            ...
                             state,            ...
                             boundary,         ...
                             mesh,             ...
                             mat,              ...
                             quadrature,       ...
                             external_source,  ...
                             fission_source )
            
           obj.d_input      = input;
           obj.d_state      = state;
           obj.d_mesh       = mesh;
           obj.d_mat        = mat;
           obj.d_quadrature = quadrature;
           obj.d_boundary   = boundary;
           
           % Setup the within-group solver.
           inner = get(input, 'inner_solver');
           if strcmp(inner, 'SourceIteration')
               obj.d_inner_solver = SourceIteration();
           elseif strcmp(inner, 'Livolant')
               obj.d_inner_solver = Livolant(); 
           end
           setup(obj.d_inner_solver,	...
                 input,                 ...
                 state,                 ...
                 boundary,              ...
                 mesh,                  ...
                 mat,                   ...
                 quadrature,            ...
                 external_source,       ...
                 fission_source);
           
        end
        
        function output = solve(obj)
            
            % Group-wise flux residual.  All residuals are *relative* and
            % are based on the L-infinity norm.
            flux_g_error = zeros(number_groups(obj.d_mat), 1);
            total_inners = 0;
            
            % Initial downscatter.  We traverse all groups, solving the
            % within-group equations.  If there is no upscatter, this
            % solves the problem.
            for g = 1:number_groups(obj.d_mat)
                
                % Set the group source.
                % group_source = d_group_s;
                
                [flux_g_error(g), inners] = ...
                    solve(obj.d_inner_solver, g);%, obj.d_source);
                
                fprintf('group %4i flux residual = %6.4e in %6i iterations \n', ...
                    g, flux_g_error(g), inners);
                
                total_inners = total_inners + inners;
            end
            
            flux_error = max(flux_g_error);
            flux_error = 1;
            iteration  = 0;
            % Outer group iterations for upscatter, if not a
            % downscatter-only problem.
            if (~downscatter(obj.d_mat))
                while (flux_error > obj.d_input.tol_outer && ...
                        iteration < obj.d_input.max_outer )
                    
                    % Iterate only on those groups for which upscattering
                    % occurs.
                    for g = 1:number_groups(obj.d_mat)
                        
                        [flux_g_error(g), inners] = ...
                            solve(obj.d_inner_solver, g);%, obj.d_source);
                        
                        fprintf('group %4i flux residual = %6.4e in %6i iterations \n', ...
                            g, flux_g_error(g), inners);
                        
                        total_inners = total_inners + inners;
                        
                    end
                    flux_error = max(flux_g_error);
                    iteration  = iteration + 1;
                end
            end
            
            output.flux_error = flux_error;
            output.flux_error = total_inners;
            
        end
        
    end
end