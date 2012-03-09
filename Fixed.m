%> @file  Fixed.m
%> @brief Fixed class definition.
% ==============================================================================
%> @brief  Fixed source multigroup solver.
%
%> The multigroup transport equations are solved via Gauss-Seidel with
%> iterations over the upscatter block. 
% ==============================================================================
classdef Fixed < handle
    
    properties
        %> User input.
        d_input
        %> Problem Geometry.
        d_mesh                 
        %> Material Database
        d_mat                   
        %> Angular mesh.
        d_quadrature            
        %> External volume source.
        d_source                
        %> External boundary source.
        d_boundary_source       
        %> Fission source.
        d_fission_source        
        %> State variables.
        d_state             
        %> Within-group solver.
        d_inner_solver 
        %> Boundary flux container.
        d_boundary
        %> Convergence tolerance.
        d_tolerance
        %> Maximum number of iterations.
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
        function this = Fixed(input,            ...
                             state,            ...
                             boundary,         ...
                             mesh,             ...
                             mat,              ...
                             quadrature,       ...
                             external_source,  ...
                             fission_source )
            
           this.d_input      = input;
           this.d_state      = state;
           this.d_mesh       = mesh;
           this.d_mat        = mat;
           this.d_quadrature = quadrature;
           this.d_boundary   = boundary;
           this.d_source     = external_source;
           this.d_fission_source = fission_source;
           
           % Check input; otherwise, set defaults for optional parameters.
           this.d_tolerance = get(input, 'outer_tolerance');
           this.d_max_iters = get(input, 'outer_max_iters');
           
           % Setup the within-group solver.
           inner = get(input, 'inner_solver');
           if strcmp(inner, 'SI')
               this.d_inner_solver = SourceIteration();
           elseif strcmp(inner, 'Livolant')
               this.d_inner_solver = Livolant(); 
           elseif strcmp(inner, 'GMRES')
               this.d_inner_solver = GMRESIteration(); 
           end
           setup(this.d_inner_solver,	...
                 input,                 ...
                 state,                 ...
                 boundary,              ...
                 mesh,                  ...
                 mat,                   ...
                 quadrature,            ...
                 external_source,       ...
                 fission_source);
           
        end
        
        function output = solve(this)
            
            % Set the boundary conditions.
            set(this.d_boundary);
            
            % Group-wise flux residual.  All residuals are *relative* and
            % are based on the L-infinity norm.
            flux_g_error = zeros(number_groups(this.d_mat), 1);
            total_inners = 0;
            
            % Initial downscatter.  We traverse all groups, solving the
            % within-group equations.  If there is no upscatter, this
            % solves the problem.
            for g = 1:number_groups(this.d_mat)
                
                [flux_g_error(g), inners] = ...
                    solve(this.d_inner_solver, g);
                
                total_inners = total_inners + inners;
                
            end
            flux_error = max(flux_g_error);
            print_iteration(this, 1, flux_error, total_inners)
            
            
            flux_error = 1;
            iteration  = 1;
            % Outer group iterations for upscatter/fission, if not a
            % downscatter-only problem.
            if (~downscatter(this.d_mat))
                
                phi_old = zeros(number_cells(this.d_mesh), ...
                                number_groups(this.d_mat));
                            
                while (flux_error > this.d_tolerance && ...
                        iteration < this.d_max_iters )
                    
                    % Save the old group fluxes
                    for g = 1:number_groups(this.d_mat)
                        phi_old(:, g) = flux(this.d_state, g);
                    end
                    
                    
                    % Iterate only on those groups for which upscattering
                    % occurs.
                    for g = upscatter_cutoff(this.d_mat):number_groups(this.d_mat)
                        
                        [flux_g_error(g), inners] = ...
                            solve(this.d_inner_solver, g);%, this.d_source);
                        
                        total_inners = total_inners + inners;
                        
                        print_iteration(this, iteration, flux_error, ...
                            total_inners)
                        
                        % Compute the norm of the group flux residual.
                        flux_g_error(g) = norm(flux(this.d_state, g) - ...
                                               phi_old(:, g));
 
                    end
                    iteration  = iteration + 1;
                    flux_error = max(flux_g_error);
                    
                end
            end
            output.flux_error   = flux_error;
            output.total_inners = total_inners;
            
        end
        
        function print_iteration(this, iteration, flux_error, total_inners)
            if (get(this.d_input, 'print_out'))
                fprintf(...
                    '-------------------------------------------------------\n')
                fprintf('       Iter: %5i, Error: %12.8f, Inners: %5i\n',...
                    iteration, flux_error, total_inners);
                fprintf(...
                    '-------------------------------------------------------\n')
            end
        end
        
        function reset_convergence(this, max_iters, tolerance)
            reset_convergence(this.d_inner_solver, max_iters, tolerance);
        end

    end
end