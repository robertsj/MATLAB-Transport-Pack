%> @file  FixedMultiply.m
%> @brief FixedMultiply class definition.
% ==============================================================================
%> @brief  Fixed source multigroup solver for multiplicative problems.
%
%> The multigroup transport equations are solved via Gauss-Seidel with
%> iterations over the upscatter block. 
% ==============================================================================
classdef FixedMultiply < Fixed
    
    properties (Access = private)
        %> User-defined keff for multiplicative problem
        d_keff     = 1;
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
        %> @return Instance of the FixedMultiply class.
        % ======================================================================
        function this = FixedMultiply(input,            ...
                                      state,            ...
                                      boundary,         ...
                                      mesh,             ...
                                      mat,              ...
                                      quadrature,       ...
                                      external_source,  ...
                                      fission_source )
                                 
            % Call base class.                                 
            this = this@Fixed(input, state, boundary, mesh, mat, ...
                              quadrature, external_source, fission_source); 
           
        end
        
        function output = solve(this)
            
            % Group-wise flux residual.  All residuals are *relative* and
            % are based on the L-infinity norm.
            flux_g_error = zeros(number_groups(this.d_mat), 1);
            total_inners = 0;
            
            % Setup the fission source.  The scaling factor (1/keff) is
            % unchanged during this iteration.  The user can set this
            % explicitly, but it defaults to unity.
            initialize(this.d_fission_source);
            setup_outer(this.d_fission_source, 1/this.d_keff);

            flux_error = 1;
            iteration  = 1;
            flag = 1;
            % Outer group iterations over all groups
            while ((flux_error > this.d_tolerance && ...
                   iteration < this.d_max_iters) || flag==1)
                
                setup_outer(this.d_fission_source, 1/this.d_keff);
               
                % Iterate over all groups
                for g = 1:number_groups(this.d_mat)
                    
                    [flux_g_error(g), inners] = ...
                        solve(this.d_inner_solver, g);
                    
                    total_inners = total_inners + inners;
                    
                    print_iteration(this, iteration, flux_error, ...
                        total_inners)
                    
                end
                
                update(this.d_fission_source);

                iteration  = iteration + 1;
                flux_error = max(flux_g_error)
                if iteration > 3
                    flag = 0;
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
        
        function set_keff(this, k)
            d_keff = keff; 
        end
    end
end