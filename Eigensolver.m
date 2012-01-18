%> @file  Eigensolver.m
%> @brief Eigensolver class definition.
% ==============================================================================
%> @brief  Solver for k-eigenvalue problems.
%
%> The standard power iteration method is employed.
% ==============================================================================
classdef Eigensolver
    
    properties
        %> User input.
        d_input
        d_mesh                  % Cartesian mesh.
        d_mat                   % Materials.
        d_quadrature            % Angular mesh.
        d_source                % External volume source.
        d_boundary_source       % External boundary source.
        d_fission_source        % Fission source.
        d_state                 % State variables.
        d_solver                % Multigroup solver.
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
        %> @return Instance of the Eigensolver class.
        % ======================================================================
        function obj = Eigensolver(  input,            ...
                                     state,            ...
                                     boundary,         ...
                                     mesh,             ...
                                     mat,              ...
                                     quadrature,       ...
                                     external_source,  ...
                                     fission_source    )
            
           obj.d_input      = input;
           obj.d_state      = state;
           obj.d_mesh       = mesh;
           obj.d_mat        = mat;
           obj.d_quadrature = quadrature;
           obj.d_boundary   = boundary;
           obj.d_fission_source = fission_source;
           
           % Setup the multigroup solver.
           obj.d_solver = ...
               Fixed(  input,            ...
                       state,            ...
                       boundary,         ...
                       mesh,             ...
                       mat,              ...
                       quadrature,       ...
                       external_source,  ...
                       fission_source );
        end
        
        function output = solve(obj)
            
            % Initialize the fission source.
            initialize(obj.d_fission_source)
            
            % Group-wise flux residual.  All residuals are *relative* and
            % are based on the L-infinity norm. "1" and "2" are "one time ago"
            % and "two times ago".
            flux_error    = 1.0; 
            flux_error_1  = 2.0;
            flux_error_2  = 3.0;
            k_eff         = 1.0;
            k_eff_1       = 1.0;
            iteration     = 0;
            
            % Do initial fission density setup
            update(obj.d_fission_source, k_eff);
            fd = density(obj.d_fission_source);
            
            % Perform outer iteration.
            while (flux_error > obj.d_input.tol_fission && ...
                   iteration < obj.d_input.max_fission )

               
                disp(' Outer iteration....')
                % Solve the multigroup equations.
                solve(obj.d_solver);
                
                % Save density and recompute.
                fd_old = fd;
                update(obj.d_fission_source, k_eff);
                fd = density(obj.d_fission_source);
                
                % Compute L-1 fission density error.
                flux_error_2 = flux_error_1;
                flux_error_1 = flux_error;
                flux_error   = max( abs(fd-fd_old)./fd );
                
                % Save and update keff
                k_eff_1 = k_eff;
                k_eff   = k_eff_1 * sum(fd) / sum(fd_old);
                
                iteration = iteration + 1;
                
                fprintf('*** Outer Iter: %5i, Error: %12.8f\n', ...
                    iteration, flux_error);
                
                if iteration > 2
                    fprintf('                  Rate: %12.8f\n', ...
                    (flux_error-flux_error_1) / (flux_error_1-flux_error_2)); 
                end
            end
            
            output.flux_error = flux_error;
            output.flux_error = total_inners;
            
        end
        
    end
end