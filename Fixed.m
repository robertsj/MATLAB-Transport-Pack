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
        d_tolerance
        d_max_iters
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
           
           % Check input; otherwise, set defaults for optional parameters.
           obj.d_tolerance = get(input, 'outer_tolerance');
           obj.d_max_iters = get(input, 'outer_max_iters');
           
           % Setup the within-group solver.
           inner = get(input, 'inner_solver');
           if strcmp(inner, 'SI')
               obj.d_inner_solver = SourceIteration();
           elseif strcmp(inner, 'Livolant')
               obj.d_inner_solver = Livolant(); 
           elseif strcmp(inner, 'GMRES')
               obj.d_inner_solver = GMRESIteration(); 
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
                
                total_inners = total_inners + inners;
            end
            flux_error = max(flux_g_error);
            print_iteration(obj, 1, flux_error, total_inners)
            
            flux_error = 1;
            iteration  = 1;
            % Outer group iterations for upscatter, if not a
            % downscatter-only problem.
            if (~downscatter(obj.d_mat))
                while (flux_error > obj.d_tolerance && ...
                        iteration < obj.d_max_iters )
                    
                    % Iterate only on those groups for which upscattering
                    % occurs.
                    for g = upscatter_cutoff(obj.d_mat):number_groups(obj.d_mat)
                        
                        [flux_g_error(g), inners] = ...
                            solve(obj.d_inner_solver, g);%, obj.d_source);
                        
                        total_inners = total_inners + inners;
                        
                        print_iteration(obj, iteration, flux_error, ...
                            total_inners)
 
                    end
                    iteration  = iteration + 1;
                    flux_error = max(flux_g_error);
                    
                end
            end
            output.flux_error   = flux_error;
            output.total_inners = total_inners;
            
        end
        
        function print_iteration(obj, iteration, flux_error, total_inners)
            if (get(obj.d_input, 'print_out'))
                fprintf(...
                    '-------------------------------------------------------\n')
                fprintf('       Iter: %5i, Error: %12.8f, Inners: %5i\n',...
                    iteration, flux_error, total_inners);
                fprintf(...
                    '-------------------------------------------------------\n')
            end
        end
        
        function reset_convergence(obj, max_iters, tolerance)
            reset_convergence(obj.d_inner_solver, max_iters, tolerance);
        end
        
    end
end