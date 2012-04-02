%> @file  InnerIteration.m
%> @brief InnerIteration class definition.
% ==============================================================================
%> @brief Solve the within-group transport equation.
%
%> The within-group transport equation in operator form is 
%> \f[
%>      \mathbf{L}\psi = 
%>      \mathbf{MS}\phi + Q
%> \f]
%> where \f$ \mathbf{L} \f$ is the streaming and collision operator,
%> \mathbf{M} is the moment-to-discrete operator, 
%> \mathbf{S} is the scattering operator, and  Q represents any 
%> source considered fixed, which includes in-scatter, fission, and
%> external sources.
%>
%> What we are really after is the scalar flux and possibly its higher
%> order moments.  Consequently, we are able to solve a somewhat different
%> problem then the within group transport equation above.  Let us operate 
%> on both sides by \f$\mathbf{L}^{-1}\f$ followed
%> by \f$ \mathbf{D}\f$ to get
%> \f[
%>      (\mathbf{I} - \mathbf{D}\mathbf{L}^{-1}\mathbf{MS})\phi
%>      = \mathbf{D} \mathbf{L}^{-1} Q \, .
%> \f] 
%> Here, \f$\mathbf{D}\f$ is the discrete-to-moment operator, defined
%> such that \f$ \phi = \mathbf{D}\psi \f$.
%>
%> Notice this is nothing but a linear system of the form 
%> \f$ \mathbf{A}x = b \f$ where
%> \f[
%>      \mathbf{A} = (\mathbf{I} - \mathbf{D}\mathbf{L}^{-1}\mathbf{MS})
%> \f] 
%> and
%> \f[
%>      b = \mathbf{D} \mathbf{L}^{-1} Q \, .
%> \f] 
%> Moreover, \f$ b\f$ is just the uncollided flux.
%>
%> A nice overview of approaches for this inner iteration is 
%> given by Larsen and Morel in <em> Nuclear Computational Science </em>.
%> We have implement the standard source iteration method, Livolant 
%> acceleration (an extrapolation technique), and solvers that use
%> MATLAB's own GMRES (and other Krylov solvers).
%> 
%> Input parameters specific to InnerIteration and derived classes:
%> - inner_max_iters (default: 100)
%> - inner_tolerance (default: 1e-5)
%>
%> \sa SourceIteration, Livolant, GMRESIteration
% ==============================================================================
classdef InnerIteration < handle
    
    properties (Access = protected)
        %> User input.
        d_input
        %> State vectors.
        d_state
        %> Boundary fluxes.
        d_boundary
        %> Problem mesh (either Cartesian mesh or MOC tracking)
        d_mesh
        %> Materials definitions.
        d_mat
        %> Angular mesh.
        d_quadrature
        %> Spatial discretization.
        d_equation
        %> Sweeper over the space-angle domain.
        d_sweeper
        %> Sweep source
        d_sweep_source
        %> User-defined external source
        d_external_source
        %> Fission source, if used
        d_fission_source
        %> Any source that remains "fixed" within the group solve
        d_fixed_source
        %> The within group scattering source
        d_scatter_source
        %> Scattering cross-section vector for faster source computation
        d_sigma_s
        %> Maximum iterations
        d_max_iters
        %> Convergence tolerance
        d_tolerance
        %
        d_g
        %> Moments to discrete operator.
        d_M
        %> Numer of sweeps
        d_number_sweeps = 0;
    end
    
    methods (Access = public)    
        
        % ======================================================================
        %> @brief Class constructor
        %
        %> @param input             Input database.
        %> @param state             State vectors, etc.
        %> @param boundary          Boundary fluxes.
        %> @param mesh              Problem mesh.
        %> @param mat               Material definitions.
        %> @param quadrature        Angular mesh..
        %> @param external_source 	User-defined external source.
        %> @param fission_source 	Fission source.
        %>
        %> @return Instance of the InnerIteration class.
        % ======================================================================
    	function this = InnerIteration()
            
        end
        
        % ======================================================================
        %> @brief Solve the within-group problem.
        %
        %> @param g     Group of the problem to be solved.
        %>
        %> @return Flux residual (L-inf norm) and iteration count.
        % ======================================================================
        [flux_error, iteration] = solve(this, g);
        
        % ======================================================================
        %> @brief Reset convergence criteria.
        %
        %> It can be useful to use a dynamically changing convergence criteria,
        %> especially for eigenproblems.  For regular power iteration, it is
        %> wasteful to use tight tolerances on the inners when the source is
        %> not as tightly known.
        %>
        %> @param max_iters     Maximum number of iterations.
        %> @param tolerance     Maximum point-wise error in flux
        % ======================================================================
        function reset_convergence(this, max_iters, tolerance)
            this.d_max_iters = max_iters;
            this.d_tolerance = tolerance;
        end
        
        function display_me(this)
            disp(' inner! ') 
        end
        
        %> @name Getters
        %> @{
        
        function n = number_sweeps(this)
            n = number_sweeps(this.d_sweeper);
        end
        
        %> @}
        
    end
    
    methods (Access = protected)
        
        % ======================================================================
        %> @brief Setup the base solver.
        %
        %> By keeping everything out of the constructor, we can avoid having to
        %> copy the signature for all derived classes, both in their
        %> definitions and whereever thisects are instantiated.  However, derived
        %> classes have other setup issues that need to be done, and so this
        %> serves as a general setup function that \em must be called before
        %> other setup stuff.
        %>
        %> @param input             Input database.
        %> @param state             State vectors, etc.
        %> @param boundary          Boundary fluxes.
        %> @param mesh              Problem mesh.
        %> @param mat               Material definitions.
        %> @param quadrature        Angular mesh..
        %> @param external_source 	User-defined external source.
        %> @param fission_source 	Fission source.
        % ======================================================================
        function this = setup_base(this,              ...
                                   input,            ...
                                   state,            ...
                                   boundary,         ...
                                   mesh,             ...
                                   mat,              ...
                                   quadrature,       ...
                                   external_source,  ...
                                   fission_source    )
                            
            this.d_input      = input;                  
            this.d_state      = state;
            this.d_boundary   = boundary;
            this.d_mat        = mat;
            this.d_mesh       = mesh;
            this.d_quadrature = quadrature;
            
            % Check input; otherwise, set defaults for optional parameters.
            this.d_tolerance = get(input, 'inner_tolerance');
            this.d_max_iters = get(input, 'inner_max_iters');
            
            % Add external source object
            this.d_external_source = external_source;
      
            % Add fission source object
            this.d_fission_source = fission_source;
            
            % Add the scattering source
            this.d_scatter_source = ScatterSource(mesh, mat, state);
            
            % Initialize the fixed source
            this.d_fixed_source = zeros(number_cells(mesh), 1);
                 
            % Initialize moment to discrete.
            this.d_M = MomentsToDiscrete(mesh.DIM);
            
            % Get discretization type.
            eq = get(input, 'equation');
            
            if mesh.DIM == 1
                
                % Equation.
                if strcmp(eq, 'DD')
                    this.d_equation = DD1D(mesh, mat, quadrature);
                elseif strcmp(eq, 'SD')
                    this.d_equation = SD1D(mesh, mat, quadrature);
                else
                   error('unsupported 1D discretization') 
                end

                this.d_sweeper = Sweep1D_mod(input, mesh, mat, quadrature, ...
                    this.d_boundary, this.d_equation);                 
                
            elseif mesh.DIM == 2 && meshed(mesh)
                
                % Equation.
                this.d_equation = DD2D(mesh, mat, quadrature);

                % Sweeper
                this.d_sweeper = Sweep2D_mod(input, mesh, mat, quadrature, ...
                    this.d_boundary, this.d_equation); 
                
            elseif mesh.DIM == 3
                
                error('3D coming soon.')
                
            elseif mesh.DIM == 2 && tracked(mesh)

                 % Equation.
                 this.d_equation = SCMOC(mesh, mat);
                 
                 % Sweeper
                 this.d_sweeper = SweepMOC(input, mesh, mat, quadrature, ...
                     this.d_boundary, this.d_equation);
                 
            else
                error('Invalid mesh dimension.')
            end
            
                  
        end

        % ======================================================================
        %> @brief Check convergence and warn if iteration limit reached.
        % ======================================================================        
        function check_convergence(this, iteration, flux_error)
            if iteration == this.d_max_iters && flux_error > this.d_tolerance
                warning('solver:convergence', ...
                    'Maximum iterations reached without convergence.')
            end
        end % end function check_convergence

        % ======================================================================
        %> @brief Print some diagnostic output.
        % ======================================================================
        function print_iteration(this, it, e0, e1, e2)
            if (get(this.d_input, 'inner_print_out'))
                fprintf('           Iter: %5i, Error: %12.8f\n', it, e0);
                if it > 2
                    fprintf('                         Rate: %12.8f\n', ...
                        (e0 - e1) / (e1 - e2));
                end
            end
        end

    end
        
end