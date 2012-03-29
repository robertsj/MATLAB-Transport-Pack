%> @file  Livolant.m
%> @brief Livolant class definition.
% ==============================================================================
%> @brief Livolant acceleration for inner iterations.
%
%> %Source iteration for a one-group transport equation is defined by the 
%> sequence
%> \f[
%>      \phi^{(l+1)} = 
%>         \frac{1}{4\pi}\mathbf{D}\mathbf{T}^{-1}\mathbf{S} \phi^{(l)}
%>       + \mathbf{D} \mathbf{T}^{-1} Q \, .
%> \f]
%> Simplifying notation, we have
%> \f[
%>      \phi^{(l+1)} = \mathbf{H} \phi^{(l)} + q\, .
%> \f]
%> %Livolant acceleration solves this via relaxation, or
%> \f[
%>      \phi^{(l+1)} = \phi^{(l)} + \mu R^{(l+1)}\, ,
%> \f]
%> where \f$ \mu \f$ is a relaxation parameter and \f$R\f$ is
%> the residual, defined as the difference between successive flux estimates.
%>
%> Here, \f$ \mu \f$ is dynamically
%> updated to minimize the residual.  In other
%> words, we choose \f$ \mu \f$ to minimize
%> \f[
%>      R^{(l+1)} = R^{(l)} + 
%>         \mu \Big ( \mathbf{H} R^{(l)}-R^{(l)} \Big ) \, .
%> \f]
%> Doing so in a least-squares sense (\em i.e. via the \f$L_2 \f$ norm) yields
%> \f[
%>     \mu^{(l)} = -\frac{ (R^{(l)})^T 
%>                        (\mathbf{H} R^{(l)}-\mathbf{R}^{(l)}) } 
%>           { || \mathbf{H} \mathbf{R}^{(l)}-\mathbf{R}^{(l)} ||^2_2 } \, .
%> \f]
%>
%> The implementation is sketched as follows (which follows Hebert).
%> Suppose we have three sequential flux estimates from regular source
%> iterations, 
%> \f[
%>      \Phi^{(0)}   =   \phi^{(l)}  \, ,
%> \f]
%> \f[
%>      \Phi^{(1)} = \mathbf{H} \Phi^{(0)} + q\, ,
%> \f]
%> and
%> \f[
%>      \Phi^{(2)} = \mathbf{H} \Phi^{(1)} + q\, .
%> \f]
%> Using these, we define the first two residuals as
%> \f[
%>      e_0 = \Phi^{(1)} - \Phi^{(0)} = R^{(l)}  \, .
%> \f]
%> \f[
%>      e_1 = \Phi^{(2)} - \Phi^{(1)} \, .
%> \f]
%> and 
%> \f[ 
%>   \mu^{(l)} = -\frac{e^T_0(e_1-e_0) } { ||e_1 - e_0)||^2_2 } \, .
%> \f]  
%> The accelerated estimate is
%> \f[
%>      \phi^{(l+1)} = \phi^{(l)} + \mu^{(l)} R^{(l)} \, .
%> \f]
%> We need three new successive estimates to reestimate \f$ \mu \f$.  These are
%> \f[
%>      \tilde{\Phi}^{(0)} = \phi^{(l+1)}  \, ,
%> \f]
%> \f[
%>      \tilde{\Phi}^{(1)} = \mathbf{H} \tilde{\Phi}^{(0)} + q =  
%>          \mu^{(l)} \Phi^{(2)} + (1-\mu^{(l)})\Phi^{(1)}\, ,
%> \f]
%> and
%> \f[
%>      \tilde{\Phi}^{(2)} = \mathbf{H} \tilde{\Phi}^{(1)} + q\, .
%> \f]
%> Notice that the first two estimates require nothing new.  Only the last
%> iterate requires further application of the operators.
%> 
%> Hebert suggests doing three source (or "free") iterations, followed
%> by three accelerated iterations for the sake of stability.  What this means
%> is that three source iterations are performed, resulting in an updated
%> \f$ \mu \f$ at a cost of three sweeps.
%> Then, the sequence of \f$\tilde{\Phi}\f$'s is performed three
%> times, with \em each sequence resulting in a new \f$ \mu \f$ at a cost of
%> a single sweep.  Hence, the Livolant iteration costs little more than source
%> iteration in computational expense and an extra flux-sized vector or two
%> in storage depending on implementation.
%>
%> The iteration limits are defaulted (to three each), but the user can change
%> these.  However, doing so may yield poor stability.  In particular,
%> negative values for \f$ \mu \f$ are problematic (and if found to be negative,
%> the accelerated iterations are simply skipped with a warning.
%>
%> Input parameters specific to Livolant:
%> - livolant_free_iters  (default: 3)
%> - livolant_accel_iters (default: 3)
%>
%> \todo Livolant, like GMRES, seems to work poorly for reflected conditions.
%>       There probably is a proper "consistent" way to include these
%>       conditions, but I don't know it yet.
%>
%> For other relevant parameters, see InnerIteration.
%>
%> Reference:
%>   Hebert, A. <em>Applied Reactor Physics</em>. See pp. 147-148 and 345-347.
%>
%> \sa SourceIteration
%>
% ==============================================================================
classdef Livolant < InnerIteration
    
    properties (Access = private)
        %> A three column vector of current, last, and penultimate fluxes.
        d_phi
        %> Column index of current flux.
        d_2 = 3;
        %> Column index of last flux.
        d_1 = 2;
        %> Column index of penultimate flux.
        d_0 = 1;
        %> Number of free iterations
        d_free_iters = 3;
        %> Number of acceleration iteration
        d_accel_iters = 3;
    end
    
    methods
       
        % ======================================================================
        %> @brief Class constructor
        %
        %> @return  Instance of the SourceIteration class.
        % ======================================================================
        function this = Livolant()
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
                                 
            % Size the working flux vector.                     
            this.d_phi = zeros(number_cells(mesh), 3);    
            
            % Set user parameters or use defaults
            if contains(input, 'livolant_free_iters')
                tmp = get(input, 'livolant_free_iters');
                if tmp > 2
                    this.d_free_iters = tmp;
                else
                    warning('user:input',...
                        ['At least 3 free iters required to update mu.', ...
                         ' Using default = 3.'])
                end
            end
            if contains(input, 'livolant_accel_iters')
                this.d_accel_iters = get(input, 'livolant_accel_iters');
            end
            
        end
        
        % ======================================================================
        %> @brief Solve the within-group problem.
        %
        %> @param g     Group of the problem to be solved.
        %>
        %> @return Flux residual (L-inf norm) and iteration (sweep) count.
        % ======================================================================
        function [flux_error, iteration] = solve(this, g)% phi, source)
            
            % Flux error currently, one time ago, and two times ago.
            flux_error    = 1.0; 
            flux_error_1  = 2.0;
            flux_error_2  = 3.0;
            iteration     = 0;
            
            % Get flux estimate for this group and put it in the current slot.
            this.d_phi(:, this.d_2) = flux(this.d_state, g);
            
            % Setup the boundary fluxes for this solve.
            set_group(this.d_boundary, g); 
            
            % Setup the equations for this group.
            setup_group(this.d_equation, g);
            
            % Build M*Q, the source that is "fixed" for this solve.
            build_fixed_source(this, g);
                 
            % Set the boundary fluxes for this sweep.
            %set(this.d_boundary);
            
            while flux_error > this.d_tolerance ...
                && iteration < this.d_max_iters
            
                % Update boundary
                update(this.d_boundary);
            
                % Perform a sequence of source ("free") iterations.
                [err, mu] = source_iterations(this, g);
                
                %fprintf('       SI Error: %12.8f,  mu: %12.8f\n', err, mu);
                iteration = iteration + this.d_free_iters;
                
                % Perform a sequence of accelerated iterations, but only if the 
                % relaxation parameter is positive (otherwise divergence lurks)
                % and we aren't already converged.
                if mu > 0 && err > this.d_tolerance
                    err = accelerated_iterations(this, g, mu);
                    %fprintf(' Livolant Error: %12.8f \n', err);
                    iteration = iteration + this.d_accel_iters;
                else
                    warning('solver:convergence', ...
                        'Negative relaxation parameter; skipping Livolant')
                end

                flux_error_2 = flux_error_1;
                flux_error_1 = flux_error;
                flux_error   = abs(err);  
                                
                print_iteration(this, iteration, flux_error, ...
                    flux_error_1, flux_error_2)

            end
           
            % Did we converge?
            check_convergence(this, iteration, flux_error);
            
            % Update the state.
            set_phi(this.d_state, this.d_phi(:, this.d_2), g);

            % Update the boundary fluxes.

        end
        
    end
    
    methods (Access = private)
        
        % ======================================================================
        %> @brief Perform a few source iterations.
        %
        %> @param   g   Current group.
        %> @param   mu  Current value of relaxation parameter.
        %> @return      Norm of current flux residual. 
        % ======================================================================
        function err = accelerated_iterations(this, g, mu)

            iteration = 0;        
            while iteration < this.d_accel_iters
                
                % First accelerated iterate.  This follows Hebert Eq. C.52.
                % We have the free iterations phi_0, phi_1, and phi_2 stored,
                % where 0 to 2 is old to new.  We overwrite the oldest flux,
                %  phi_0 = mu*phi_1 + (1-mu)*phi_0.
                this.d_phi(:, this.d_0) = mu*this.d_phi(:, this.d_1) + ...
                    (1 - mu)*this.d_phi(:, this.d_0);
                
                % Second accelerated iterate.  This follows Hebert Eq. C.53.  
                % Again, we need only previous fluxes and mu.
                this.d_phi(:, this.d_1) = mu*this.d_phi(:, this.d_2) + ...
                    (1 - mu)*this.d_phi(:, this.d_1);
                
                % Compute the first residual.
                e0 = this.d_phi(:, this.d_1) - this.d_phi(:, this.d_0);
                
                % **Now** we need to do a sweep for the third flux.
                
                % Build the within-group source (using phi_1!)
                build_scatter_source(this, g, this.d_phi(:, this.d_1));
                
                % Construct total sweep source.
                sweep_source = (this.d_fixed_source + this.d_scatter_source); 
                
                % Set incident boundary fluxes.
                %set(this.d_boundary);
                
                % Sweep over all angles and meshes.            
                this.d_phi(:, this.d_2) = sweep(this.d_sweeper, sweep_source, g);      
                
                % Compute the second residual
                e1 = this.d_phi(:, this.d_2) - this.d_phi(:, this.d_1);
                
                % err = max( abs(e1)./this.d_phi(:, this.d_2) );
                err = norm(e1);
                
                % Update the relaxation factor.
                mu = relaxation_factor(this, e0, e1);
                
                iteration = iteration + 1;
                   
            end
            
        end % end function accelerated_iterations
        
        % ======================================================================
        %> @brief Perform a few source iterations.
        %
        %> @param   g   Current group.
        %>
        %> @return      Relaxation parameter and norm of flux residual. 
        % ======================================================================
        function [err, mu] = source_iterations(this, g)

            iteration = 0; 
            while iteration < this.d_free_iters
            
                % Build the within-group source
                build_scatter_source(this, g, this.d_phi(:, this.d_2));
                
                % Construct total sweep source.
                sweep_source = (this.d_fixed_source + this.d_scatter_source); 
                
                % Set incident boundary fluxes.
                % set(this.d_boundary);
                
                % Shuffle indices so new flux is put in right storage.
                shuffle_index(this);
                
                % Sweep over all angles and meshes.            
                this.d_phi(:, this.d_2) = sweep(this.d_sweeper, sweep_source, g);            
                
                iteration = iteration + 1;   
                % Compute residuals following the last two iterations
                if iteration == this.d_free_iters - 1
                    e0 = this.d_phi(:, this.d_2) - this.d_phi(:, this.d_1);
                elseif iteration == this.d_free_iters 
                    e1 = this.d_phi(:, this.d_2) - this.d_phi(:, this.d_1);
                    % err = max( abs(e1)./this.d_phi(:, this.d_2) );
                    err = norm(e1);
                end
                      
            end

            % Compute relaxation factor
            mu = relaxation_factor(this, e0, e1);
            
        end % end function source_iterations
        
        % ======================================================================
        %> @brief Shuffle indices so that storage is reused.
        %
        %> This is performed \em before and iteration so that the "2" flux can
        %> be written into.
        % ======================================================================
        function this = shuffle_index(this)
            temp    = this.d_0;
            this.d_0 = this.d_1;
            this.d_1 = this.d_2;
            this.d_2 = temp;
        end % end function shuffle_index
        
        % ======================================================================
        %> @brief Compute the relaxation factor from two residual vectors.
        %
        %> @param   e0      Previous residual.
        %> @param   e1      Current residual.
        %> @return          Relaxation factor, \f$ \mu \f$.
        % ======================================================================
        function mu = relaxation_factor(this, e0, e1)
            del_e       = e1 - e0;
            del_e_sq    = del_e'*del_e;
            mu          = -(e0' * del_e) / del_e_sq;  
        end % end function relaxation_factor
        
    end
    
end