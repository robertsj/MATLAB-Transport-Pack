%> @file  Eigensolver.m
%> @brief Eigensolver class definition.
% ==============================================================================
%> @brief  Solver for k-eigenvalue problems.
%
%> The homogeneous transport equation in operator form is
%> \f[
%>      \mathbf{T}\psi = 
%>        \mathbf{MSD}\psi + \frac{1}{k} \mathbf{MFD} \psi \, ,
%> \f]
%>
%> 
%> The standard power iteration method is employed.
% ==============================================================================
classdef Eigensolver < handle
    
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
        %
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
        %> @return Instance of the Eigensolver class.
        % ======================================================================
        function this = Eigensolver(  input,            ...
                                     state,            ...
                                     boundary,         ...
                                     mesh,             ...
                                     mat,              ...
                                     quadrature,       ...
                                     external_source,  ...
                                     fission_source    )
            
           this.d_input      = input;
           this.d_state      = state;
           this.d_mesh       = mesh;
           this.d_mat        = mat;
           this.d_quadrature = quadrature;
           this.d_boundary   = boundary;
           this.d_fission_source = fission_source;
           
           % Check input; otherwise, set defaults for optional parameters.
           this.d_tolerance = get(input, 'eigen_tolerance');
           this.d_max_iters = get(input, 'eigen_max_iters');
           
           % Setup the multigroup solver.
           this.d_solver = ...
               Fixed(  input,            ...
                       state,            ...
                       boundary,         ...
                       mesh,             ...
                       mat,              ...
                       quadrature,       ...
                       external_source,  ...
                       fission_source );
        end
        
        function output = solve(this)
            
            disp('*** Beginning eigensolve***')
            
            % Initialize the fission source.
            initialize(this.d_fission_source);
            
            % Initialize the flux to unity in all groups.
            phi = ones(number_cells(this.d_mesh), 1);
            for g = 1:number_groups(this.d_mat)
               set_phi(this.d_state, phi, g);
            end
            
            % Do initial fission density setup
            update(this.d_fission_source);
            
            % Normalize the density to unity
            normalize(this.d_fission_source);
            fd = density(this.d_fission_source);
            
            % Group-wise flux residual.  All residuals are *relative* and
            % are based on the L-infinity norm. "1" and "2" are "one time ago"
            % and "two times ago".
            flux_error    = 1.0; 
            flux_error_1  = 2.0;
            flux_error_2  = 3.0;
            k_eff         = 1.0;
            k_eff_1       = 1.0;
            iteration     = 0;

            sweeps = 0;
            
            % Perform outer iteration.
            while (flux_error > this.d_tolerance && ...
                   iteration  < this.d_max_iters )
               
                % Setup the fission source.
                setup_outer(this.d_fission_source, 1/k_eff);
               
                % Solve the multigroup equations.
                out     = solve(this.d_solver);
                
                % Track number of sweeps.
                sweeps  = sweeps + out.total_inners;
                
                % Save density and recompute.
                update(this.d_fission_source);
                fd_old = fd;
                fd     = density(this.d_fission_source);
                
                % Compute L-1 fission density error.
                flux_error_2 = flux_error_1;
                flux_error_1 = flux_error;
                flux_error   = max( abs(fd-fd_old)./fd );
                
                %reset_convergence(this.d_solver, 100, flux_error/10);
                
                % Save and update keff
                k_eff_1 = k_eff;
                k_eff   = k_eff_1 * max(fd) / max(fd_old);
                
                iteration = iteration + 1;
                print_iteration(this, iteration, sweeps, ...
                    flux_error, flux_error_1, flux_error_2, k_eff)
            end
            set_eigenvalue(this.d_state, k_eff);
            output.flux_error = flux_error;
            output.sweeps = sweeps;
            
        end
        
        
        function print_iteration(this, it, sweeps, e0, e1, e2, k)
            fprintf('=========================================================================\n')
            fprintf('    Iter: %5i, Error: %12.8f, keff: %10.6f, Sweeps: %5i\n',...
                it, e0, k, sweeps);
            if it > 2
                fprintf('                      Rate: %12.8f\n', ...
                    (e0 - e1) / (e1 - e2));
            end
            fprintf('=========================================================================\n')
        end
        
    end
end