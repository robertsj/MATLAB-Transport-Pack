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
        %>
        d_input
        d_state
        d_boundary
        d_mesh
        d_mat
        d_quadrature
        d_equation
        d_sweeper
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
        % 
        d_max_iters
        d_tolerance
        %
        d_g
        %> Moments to discrete operator.
        d_M
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
    	function obj = InnerIteration()
            % Nothing here for now.
        end
        
        % ======================================================================
        %> @brief Solve the within-group problem.
        %
        %> @param g     Group of the problem to be solved.
        %>
        %> @return Flux residual (L-inf norm) and iteration count.
        % ======================================================================
        [flux_error, iteration] = solve(obj, g);
        
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
        function reset_convergence(obj, max_iters, tolerance)
            obj.d_max_iters = max_iters;
            obj.d_tolerance = tolerance;
        end
    end
    
    methods (Access = protected)
        
        % ======================================================================
        %> @brief Setup the base solver.
        %
        %> By keeping everything out of the constructor, we can avoid having to
        %> copy the signature for all derived classes, both in their
        %> definitions and whereever objects are instantiated.  However, derived
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
        function obj = setup_base(obj,              ...
                                  input,            ...
                                  state,            ...
                                  boundary,         ...
                                  mesh,             ...
                                  mat,              ...
                                  quadrature,       ...
                                  external_source,  ...
                                  fission_source    )
                            
            obj.d_state      = state;
            obj.d_boundary   = boundary;
            obj.d_mat        = mat;
            obj.d_mesh       = mesh;
            obj.d_quadrature = quadrature;
            
            % Check input; otherwise, set defaults for optional parameters.
            obj.d_tolerance = get(input, 'inner_tolerance');
            obj.d_max_iters = get(input, 'inner_max_iters');
            
            % Add external source
            obj.d_external_source = external_source;
      
            % Add fission source
            obj.d_fission_source = fission_source;
              
            % Initialize the group source
            obj.d_fixed_source = zeros(number_cells(mesh), 1);
            
            % Initialize the within-group source
            obj.d_scatter_source = zeros(number_cells(mesh), 1);
            
            % Initialize scatter data.
            initialize_scatter(obj);
            
            % Initialize moment to discrete.
            obj.d_M = MomentsToDiscrete(mesh.DIM);
            
            if mesh.DIM == 1
                
                % Equation.
                obj.d_equation = DD1D(mesh, mat);

                % Sweeper
                obj.d_sweeper = Sweep1D(input, mesh, mat, quadrature, ...
                    obj.d_boundary, obj.d_equation); 
                
            elseif mesh.DIM == 2
                
                % Equation.
                obj.d_equation = DD2D(mesh, mat);

                % Sweeper
                obj.d_sweeper = Sweep2D(input, mesh, mat, quadrature, ...
                    obj.d_boundary, obj.d_equation); 
                
            elseif mesh.DIM == 3
                
            elseif mesh.MOC 
                 % Example of how MOC can be plugges right in.
            else
                error('Invalid mesh dimension.')
            end
            
                  
        end
        
        
        % ======================================================================
        %> @brief Prebuild scattering matrix for each cell.
        %
        %> While not the most *memory* efficient, this saves *time* by
        %> eliminate lots of loops.
        % ======================================================================
        function obj = initialize_scatter(obj)

            for g = 1:number_groups(obj.d_mat)
            	obj.d_sigma_s{g} = zeros(number_cells(obj.d_mesh), ...
                                         number_groups(obj.d_mat));
            end
            
            mat = reshape(mesh_map(obj.d_mesh, 'MATERIAL'), ...
                number_cells(obj.d_mesh), 1);
            
            for i = 1:number_cells(obj.d_mesh)
                for g = 1:number_groups(obj.d_mat)
                    for gp = lower(obj.d_mat, g):upper(obj.d_mat, g)
                        obj.d_sigma_s{g}(i, gp) = ...
                            sigma_s(obj.d_mat, mat(i), g, gp);
                    end
                end
            end
  
        end % end function initialize_scatter
        
        % ======================================================================
        %> @brief Build the within-group scattering source.
        %
        %> @param   g       Group for this problem.
        %> @param   phi     Current group flux.
        % ======================================================================
        function obj = build_scatter_source(obj, g, phi)

            % Build the source.
            obj.d_scatter_source = ...
                phi .* obj.d_sigma_s{g}(:, g);
            
            % Apply moments-to-discrete operator.
            obj.d_scatter_source = apply(obj.d_M, obj.d_scatter_source); 
            
        end % end function build_scatter_source
        
        % ======================================================================
        %> @brief Build source for this group excluding within-group scatter.
        %
        %> @param   g       Group for this problem.
        % ======================================================================
        function obj = build_fixed_source(obj, g)
        % function Q = build_source(obj, g)
        %   This builds the source (excluding within-group scatter) for the
        %   sweep.  Specifically, a *discrete* source is generated, i.e. the
        %   source appropriate for T*psi = Q.
            obj.d_g = g;
            fprintf('          Group: %5i\n', g);
            
            obj.d_fixed_source(:) = 0;
            
            % Add downscatter source.
            for gp = lower(obj.d_mat, g) : g - 1
                
                % Get the group gp flux.
                phi = flux(obj.d_state, gp);

                % Add group contribution.
                obj.d_fixed_source = obj.d_fixed_source + ...
                	 phi .* obj.d_sigma_s{g}(:, gp);
                 
            end
                
            % Add upscatter source. 
            for gp = g + 1 : upper(obj.d_mat, g)
                
                % Get the group gp flux.
                phi = flux(obj.d_state, gp);
                

                % Add group contribution.
                obj.d_fixed_source = obj.d_fixed_source + ...
                	 phi .* obj.d_sigma_s{g}(:, gp);
                 
            end
            
            % Apply the moments-to-discrete operator.
            obj.d_fixed_source = apply(obj.d_M, obj.d_fixed_source);
            
            % Add the fission source if required.
            if (initialized(obj.d_fission_source))

                % Get the group gp fission source.
                f = source(obj.d_fission_source, g);
                
                % Add it.  This *assumes* the fission source returns a
                % vector prescaled to serve as a discrete source.
                obj.d_fixed_source = obj.d_fixed_source + f;
                
            end

            
            % External (if required)
            if (initialized(obj.d_external_source))
                
            	obj.d_fixed_source = obj.d_fixed_source + ...
                    source(obj.d_external_source, g);
                
            end
            
        end % end function build_fixed_source
        
        function check_convergence(obj, iteration, flux_error)
            if iteration == obj.d_max_iters && flux_error > obj.d_tolerance
                warning('solver:convergence', ...
                    'Maximum iterations reached without convergence.')
            end
        end % end function check_convergence

    
        function print_iteration(obj, it, e0, e1, e2)
            fprintf('           Iter: %5i, Error: %12.8f\n', it, e0);
            if it > 2
                fprintf('                         Rate: %12.8f\n', ...
                    (e0 - e1) / (e1 - e2));
            end
        end
       
        
    end
        
end