%> @file  GMRESIteration.m
%> @brief GMRESIteration class definition.
% ==============================================================================
%> @brief %Source iteration for inner iterations.
%
%> %Source iteration for a one-group transport equation is defined by the 
%> sequence
%> \f[
%>      \phi^{(l+1)} = 
%>         \frac{1}{4\pi}\mathbf{D}\mathbf{T}^{-1}\mathbf{S} \phi^{(l)}
%>       + \mathbf{D} \mathbf{T}^{-1} Q \, .
%> \f]
%>
%> In many cases this algorithm is too slow, particularly for high scattering, 
%> highly diffusive problems.  In these cases, an acceleration scheme is
%> required; these schemes include the simple \ref Livolant and \ref 
%> GMRES methods that work directly with the system above.  
%> Ultimately, these are not 
%> "acceleration" schemes in the standard sense of the word but rather are
%> different (and more efficient) approaches to solving the linear system 
%> defined by the one group
%> transport equation.  
%>
%> Other acceleration schemes that work within the source iteration include 
%> diffusion synthetic acceleration (DSA), coarse mesh finite difference (CMFD)
%> acceleration, and coarse mesh rebalance (CMR).  Each of these will be
%> implemented in time.
% ==============================================================================
classdef GMRESIteration < InnerIteration
    
    methods
       
        % ======================================================================
        %> @brief Class constructor
        %
        %> @return  Instance of the GMRESIteration class.
        % ======================================================================
        function obj = GMRESIteration()
            % Nothing here for now.
        end
        
        % ======================================================================
        %> @brief Setup the solver.
        %
        %> @param input             Input database.
        %> @param state             State vectors, etc.
        %> @param boundary          Boundary fluxes.
        %> @param mesh              Problem mesh.
        %> @param mat               Material definitions.
        %> @param quadrature        Angular mesh..
        %> @param external_source 	User-defined external source.
        %> @param fission_source 	Fission source.
        % ======================================================================
        function obj = setup(obj,              ...
                             input,            ...
                             state,            ...
                             boundary,         ...
                             mesh,             ...
                             mat,              ...
                             quadrature,       ...
                             external_source,  ...
                             fission_source    )
                         
            % First do base class setup.
            setup_base( obj,              ...
                        input,            ...
                        state,            ...
                        boundary,         ...
                        mesh,             ...
                        mat,              ...
                        quadrature,       ...
                        external_source,  ...
                        fission_source);
            
            % Nothing else here for now.
        end
        
        % ======================================================================
        %> @brief Solve the within-group problem.
        %
        %> @param g     Group of the problem to be solved.
        %>
        %> @return Flux residual (L-inf norm) and iteration count.
        % ======================================================================
        function [flux_error, iteration] = solve(obj, g)% phi, source)

            % Setup the boundary fluxes for this solve.
            initialize(obj.d_boundary, g);
            
            % Setup the equations for this group.
            setup_group(obj.d_equation, g);
            
            % Build the fixed source.
            build_fixed_source(obj, g);
            
            % Compute the uncollided flux, i.e. phi_uc = D*inv(T)*Q_fixed.
            % This is the right hand side.
            B = sweep(obj.d_sweeper, obj.d_fixed_source, g); 
            

            % Set the initial flux
            %phi = flux(obj.d_state, g);
            
            % Call GMRES
            [phi, flag, flux_error, iter] = ...
                gmres(@(x)apply(x, obj), B, 14, obj.d_tolerance, ...
                obj.d_max_iters, [], [], B);
            fprintf('           GMRES Outers: %5i,  Inners: %5i\n', ...
                iter(1), iter(2));
            switch flag
                case 0
                    % Okay.
                case 1
                     warning('solver:convergence', ...
                         'GMRES iterated MAXIT times without converging.')
                case 2
                     warning('solver:convergence', ...
                         'GMRES preconditioner was ill-conditioned.')  
                case 3
                     warning('solver:convergence', ...
                         'GMRES stagnated.')
                otherwise
                    error('GMRES returned unknown flag.')
            end
            
            iteration = iter(1) * iter(2);
            
            % Did we converge?
            check_convergence(obj, iteration, flux_error);
            
            % Update the state.
            set_phi(obj.d_state, phi, g);

            % Update the boundary fluxes.

        end
        
    end
end

function y = apply(x, obj)

    % Build the within-group source.  This is equivalent to
    %   x <-- M*S*x
    build_scatter_source(obj, obj.d_g, x); 
    sweep_source = obj.d_scatter_source;
    
    % Set incident boundary fluxes.
    set(obj.d_boundary);
    
    % Sweep over all angles and meshes.  This is equivalent to
    %   y <-- D*inv(L)*M*S*x
    y = sweep(obj.d_sweeper, sweep_source, obj.d_g);
    
    % Now, return the following
    % y <-- x - D*inv(L)*M*S*x = (I - D*inv(L)*M*S)*x 
    y = x - y;
    
end