%> @file  SourceIteration.m
%> @brief SourceIteration class definition.
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
classdef SourceIteration < InnerIteration
    
    methods
       
        % ======================================================================
        %> @brief Class constructor
        %
        %> @return  Instance of the SourceIteration class.
        % ======================================================================
        function this = SourceIteration()
            
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
            
            % Flux error currently, one time ago, and two times ago.
            flux_error    = 1.0; 
            flux_error_1  = 2.0;
            flux_error_2  = 3.0;
            iteration  = 0;
              
            % Setup the boundary fluxes for this solve.
            set_group(this.d_boundary, g);
            
            % Setup the equations for this group.
            setup_group(this.d_equation, g);
            
            % Build the fixed source.  Since it's unchanging, do it once.
            build_fixed_source(this, g);

            % Set the initial flux
            phi      = flux(this.d_state, g);
            phi_old  = phi;
            
            while flux_error > this.d_tolerance && iteration < this.d_max_iters
            
                % Construct total sweep source.  This is just the sum of the
                % within-group and other terms.  These are *moments*, and
                % so the sweeper has to get the discrete component for an
                % angle if needed.
                sweep_source = build_sweep_source(this, g, phi);
                
                % Update incident boundary fluxes.
                update(this.d_boundary);

                % Sweep over all angles and meshes.
                phi = sweep(this.d_sweeper, sweep_source, g);            

                % Convergence criteria and diagnostics
                flux_error_2 = flux_error_1;
                flux_error_1 = flux_error;

                % Relative difference might be better.
                flux_error   = norm(phi-phi_old);  
                
                phi_old      = phi;
                iteration    = iteration + 1;   

                if (mod(iteration, 20)==0)
                    print_iteration(this, iteration, flux_error, ...
                        flux_error_1, flux_error_2)
                end
                
            end
            
            % Did we converge?
            check_convergence(this, iteration, flux_error);
            
            % Update the state.
            set_phi(this.d_state, phi, g);

        end
        
    end
    
    methods (Access = protected)
       
        % ======================================================================
        %> @brief Build the fixed source.
        % ======================================================================
        function this = build_fixed_source(this, g)
            
            % Add in-scatter source.
            this.d_fixed_source = ...
                build_in_scatter_source(this.d_scatter_source, g);
        	
            % Add the fission source if present.
            if (initialized(this.d_fission_source))
                this.d_fixed_source = this.d_fixed_source + ...
                                      source(this.d_fission_source, g);
            end
     
            % Add the external source if present.
            if (initialized(this.d_external_source))
            	this.d_fixed_source = this.d_fixed_source + ...
                                      source(this.d_external_source, g);   
            end
            
        end  
        
        % ======================================================================
        %> @brief Build and return the sweep source.
        % ======================================================================
        function q = build_sweep_source(this, g, phi)
           
        	% Build the within-group source, S_gg * phi_g
            q = build_within_group_source(this.d_scatter_source, g, phi);
            
            % Add fixed source component.
            q = q + this.d_fixed_source;
            
            % Apply moments-to-discrete operator
            %q = apply(this.d_M, q);
            
        end
 
    end
    
end