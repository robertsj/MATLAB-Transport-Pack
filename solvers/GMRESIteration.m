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
%>
%> @todo Currently, GMRES doesn't like reflecting boundaries.  I'm sure
%>       I'm just missing something...
% ==============================================================================
classdef GMRESIteration < InnerIteration
    
    properties
        %> Krylov solver type.
        d_krylov 
        %>
        d_diffop
        %>
        d_pc = 0
        d_apply_m = []
    end
    
    methods
       
        % ======================================================================
        %> @brief Class constructor
        %
        %> @return  Instance of the GMRESIteration class.
        % ======================================================================
        function this = GMRESIteration()
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
        function this = setup(this,              ...
                             input,            ...
                             state,            ...
                             boundary,         ...
                             mesh,             ...
                             mat,              ...
                             quadrature,       ...
                             external_source,  ...
                             fission_source    )
                         
            % First do base class setup.
            setup_base( this,              ...
                        input,            ...
                        state,            ...
                        boundary,         ...
                        mesh,             ...
                        mat,              ...
                        quadrature,       ...
                        external_source,  ...
                        fission_source);
                    
            % Set user parameters or use defaults
            if input.get('inner_krylov_solver')
                this.d_krylov = input.get('inner_krylov_solver');
            else
                this.d_krylov = 'gmres';
            end
            if input.get('inner_precondition')
                this.d_pc = 1;
                this.d_diffop = ...
                    DiffusionOperator(input, mat, mesh);
                this.d_apply_m = @(x)apply_m(x, this);
            end
            % Nothing else here for now.
        end
        
        % ======================================================================
        %> @brief Solve the within-group problem.
        %
        %> @param g     Group of the problem to be solved.
        %>
        %> @return Flux residual (L-inf norm) and iteration count.
        % ======================================================================
        function [flux_error, iteration] = solve(this, g)% phi, source)

            % Setup the boundary fluxes for this solve.
            set(this.d_boundary);          % Update the conditions.
            set_group(this.d_boundary, g); % Set the group

            % Setup the equations for this group.
            setup_group(this.d_equation, g);
            
            % Build the fixed source.
            build_fixed_source(this, g);
            
            % B = D*inv(L)*Q
            B = sweep(this.d_sweeper, this.d_fixed_source, g); 
            
            % Unset the boundary.
            reset(this.d_boundary);       % Conditions now homogeneous
            
            print_out = get(this.d_input, 'inner_print_out');
            
            % Call the Krylov solver
            if strcmp(this.d_krylov, 'gmres')
                [phi, flag, flux_error, iter] = ...
                    gmres(@(x)apply(x, this), B, 30, this.d_tolerance, ...
                    40, this.d_apply_m, [], B);
                if print_out
                    fprintf('           GMRES Outers: %5i,  Inners: %5i\n', ...
                        iter(1), iter(2));
                end
                iteration = iter(1) * iter(2);
            elseif strcmp(this.d_krylov, 'bicgstab')
                [phi, flag, flux_error, iter] = ...
                    bicgstab(@(x)apply(x, this), B, this.d_tolerance, ...
                    40, this.d_apply_m, [], B);
                if print_out
                    fprintf('         BiCGStab Iters: %5i \n', iter);
                end
                iteration = iter;
            elseif strcmp(this.d_krylov, 'bicgstabl')
                [phi, flag, flux_error, iter] = ...
                    bicgstabl(@(x)apply(x, this), B, this.d_tolerance, ...
                    40, this.d_apply_m, [], B);
                if print_out
                    fprintf('      BiCGStab(l) Iters: %5i \n', iter);
                end
                iteration = iter;
            else
                error('Invalid Krylov solver selected for inners.')
            end
            
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

            % Did we converge?
            check_convergence(this, iteration, flux_error);
            
            % Update the state.
            set_phi(this.d_state, phi, g);

            % Update the boundary fluxes.

        end
        
    end
end

function y = apply(x, this)

    % Build the within-group source.  This is equivalent to
    %   x <-- M*S*x
    build_scatter_source(this, this.d_g, x); 
        
    sweep_source = this.d_scatter_source;
    
    % Set incident boundary fluxes.  Note, the set command loops through all
    % groups, and hence we must explicitly set the group after.
    set(this.d_boundary);
    set_group(this.d_boundary, this.d_g); 
    
    % Sweep over all angles and meshes.  This is equivalent to
    %   y <-- D*inv(L)*M*S*x
    y = sweep(this.d_sweeper, sweep_source, this.d_g);  
    
    % Now, return the following
    % y <-- x - D*inv(L)*M*S*x = (I - D*inv(L)*M*S)*x 
    y = x - y;

end

%> @brief Apply one-group diffusion preconditioner.
function y = apply_m(x, this)

    M = this.d_diffop.get_1g_operator(this.d_g);
    % y <-  (I + inv(C)*S)x 
    build_scatter_source(this, this.d_g, x, 0); % undo mom-to-dis
    y = x + M\this.d_scatter_source;

end