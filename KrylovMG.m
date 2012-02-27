%> @file  KrylovMG.m
%> @brief KrylovMG class definition.
% ==============================================================================
%> @brief  Krylov multigroup solver.
%
%> The multigroup transport equations are solved via a Krylov method.  The
%> user can specify whether the fission source is to be updated (as would
%> desirable for fixed source problems with multiplication)
% ==============================================================================
classdef KrylovMG < InnerIteration
    
    properties
%         %> User input.
%         d_input
%         %> Problem Geometry.
%         d_mesh                 
%         %> Material Database
%         d_mat                   
%         %> Angular mesh.
%         d_quadrature            
%         %> External volume source.
%         d_source                
%         %> External boundary source.
%         d_boundary_source       
%         %> Fission source.
%         d_fission_source        
%         %> State variables.
%         d_state             
%         %> Within-group solver.
%         d_inner_solver 
%         %> Boundary flux container.
%         d_boundary
%         %> Convergence tolerance.
%         d_tolerance
%         %> Maximum number of iterations.
%         d_max_iters
    end
    
    methods
        
        % ======================================================================
        %> @brief Class constructor
        %
        %> @param input             User input.
        %> @param state             State vectors.
        %> @param boundary          Boundary flux container.
        %> @param mesh              Geometry.
        %> @param mat               Material database.
        %> @param quadrature        Angular mesh.
        %> @param external_source 	Fixed source.
        %> @param fission_source 	Fission source.
        %>
        %> @return Instance of the KrylovMG class.
        % ======================================================================
        function this = KrylovMG(input,            ...
                                state,            ...
                                boundary,         ...
                                mesh,             ...
                                mat,              ...
                                quadrature,       ...
                                external_source,  ...
                                fission_source )
            
            % First do base class setup.  This builds the scattering
            % matrices, etc.
            setup_base( this,              ...
                        input,            ...
                        state,            ...
                        boundary,         ...
                        mesh,             ...
                        mat,              ...
                        quadrature,       ...
                        external_source,  ...
                        fission_source);
           
        end
        
        % ======================================================================
        %> @brief Solve the multigroup fixed source problem.
        %> @return Instance of the KrylovMG class.
        % ======================================================================
        function output = solve(this)
            
            % Setup.
            n = number_cells(this.d_mesh);
            for g = 1:number_groups(this.d_mat)
                % Define boundary fluxes.
                initialize(this.d_boundary, g);
                % Setup the equations for this group.
                setup_group(this.d_equation, g);
                % Build the fixed source.
                build_fixed_source(this, g);
                % Compute the uncollided flux (i.e. RHS)
                B((g-1)*n+1:g*n, 1) = ...
                    sweep(this.d_sweeper, this.d_fixed_source, g); 
            end
            
            [phi, flag, flux_error, iter] = ...
                gmres(@(x)apply(x, this),   ... % Function to apply operator
                B,                          ... % right hand side
                30,                         ... % restart
                this.d_tolerance,           ... % tolerance
                40,                         ... % maxit (maxit*restart = total # applications)
                [],                         ... % left pc
                [],                         ... % right pc
                B                           );  % initial guess
            
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
            
            fprintf('        MG GMRES Outers: %5i,  Inners: %5i\n', ...
                iter(1), iter(2));
            iteration = iter(1) * iter(2);
            
            output.flux_error   = flux_error;
            output.total_inners = iteration;
        end
    end
    
    methods (Access = protected)
        
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

    end
end


% ======================================================================
%> @brief Apply the multigroup transport operator.
%> @return Matrix-vector
% ======================================================================
function y = apply(x, this)


    % Number of unknowns per group
    n   = number_cells(this.d_mesh);
    ng  = number_groups(this.d_mat);
    % Store the incoming Krylov vector
    phi = reshape(x,n,ng);
    y   = phi;
    % Build the application of 
    %
    %   |   I-T*M*S_11    -T*M*S_12    -T*M*S_13   ... | |phi_1|   |b_1|
    %   |    -T*M*S_21   I-T*M*S_22    -T*M*S_23   ... |*|phi_2| = |b_2|
    %   |    ...                                       | |phi_3|   |b_3|
    %  
    % where T = D*inv(L), D is the discrete to moments operator, and inv(L)
    % is the space-angle sweep.
    for g = 1:ng

        % Build all scattering sources.
        build_all_scatter_source(this, g, phi);
        
        % This is the *only* contribution to this groups's sweep source
        sweep_source = this.d_scatter_source;

        % Set incident boundary fluxes.
        set(this.d_boundary);

        % Sweep over all angles and meshes.  This is equivalent to
        %   y <-- D*inv(L)*M*S*x
        y(:, g) = sweep(this.d_sweeper, sweep_source, this.d_g);

        % Now, return the following
        % y <-- x - D*inv(L)*M*S*x = (I - D*inv(L)*M*S)*x 
        y(:, g) = phi(:, g) - y(:, g);
        
    end
    y = reshape(y, n*ng, 1);
    
end