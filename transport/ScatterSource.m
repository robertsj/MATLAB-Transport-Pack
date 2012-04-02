%> @file  ScatterSource.m
%> @brief ScatterSource class definition.
% ==============================================================================
%> @brief Compute the scattering source for various transport problems.
%
%> 
% ==============================================================================
classdef ScatterSource < handle
    
    properties (Access = protected)
        %> Mesh
        d_mesh
        %> Materials
        d_mat
        %> State vectors
        d_state
        %> Scattering cross section vector (for fast computation)
        d_sigma_s
    end
    
    methods (Access = public)
       
        function this = ScatterSource(mesh, mat, state)
            this.d_mesh     = mesh;
            this.d_mat      = mat;
            this.d_state    = state;
            % Initialize
            this.initialize_scatter();
        end
       
        % ======================================================================
        %> @brief Prebuild scattering matrix for each cell.
        %
        %> While not the most *memory* efficient, this saves *time* by
        %> eliminate lots of loops.
        %>
        %> This should be called just once *before* any solves.
        % ======================================================================
        function this = initialize_scatter(this)

            for g = 1:number_groups(this.d_mat)
            	this.d_sigma_s{g} = zeros(number_cells(this.d_mesh), ...
                                         number_groups(this.d_mat));
            end
            
            % Get the fine mesh material map or the region material map.
            if meshed(this.d_mesh)
                mat = reshape(mesh_map(this.d_mesh, 'MATERIAL'), ...
                    number_cells(this.d_mesh), 1);
            else                
                mat = region_mat_map(this.d_mesh);
            end
            
            % Build the local scattering matrix.
            for i = 1:number_cells(this.d_mesh)
                for g = 1:number_groups(this.d_mat)
                    for gp = lower(this.d_mat, g):upper(this.d_mat, g)
                        this.d_sigma_s{g}(i, gp) = ...
                            sigma_s(this.d_mat, mat(i), g, gp);
                    end
                end
            end

        end % end function initialize_scatter
        
        % ======================================================================
        %> @brief Build the within-group scattering source.
        %>
        %> This constructs
        %> @f[
        %>     q_g = \mathbf{S}_{gg}\phi_g \, .
        %> @f]
        %>
        %> @param   g       Group for this problem.
        %> @param   phi     Group flux.
        %> @return          Within-group catter source.
        % ======================================================================
        function q = build_within_group_source(this, g, phi)
            q = phi .* this.d_sigma_s{g}(:, g);
        end % end function build_scatter_source
        
        % ======================================================================
        %> @brief Build the in-scatter source.
        %>
        %> This is likely used to construct the inscatter source for a
        %  within-group solve, for which inscatter is viewed as "fixed".
        %> This constructs
        %> @f[
        %>     q_g = \sum^G_{g',g\ne g'} \mathbf{S}_{gg'}\phi_{g'} \, .
        %> @f]
        %>
        %> This assumes the state is up-to-date.
        %> 
        %> @param   g       Group for this problem.
        %> @return          Within-group catter source.
        % ======================================================================
        function q = build_in_scatter_source(this, g)
            
            % Add downscatter source.
            for gp = lower(this.d_mat, g) : g - 1        
                % Get the group gp flux.
                phi = flux(this.d_state, gp);
                % Add group contribution.
                q = q + phi .* this.d_sigma_s{g}(:, gp);
            end
                
            % Add upscatter source. 
            for gp = g + 1 : upper(this.d_mat, g)
                % Get the group gp flux.
                phi = flux(this.d_state, gp);
                % Add group contribution.
                q = q + phi .* this.d_sigma_s{g}(:, gp);
            end
            
        end % end function build_scatter_source
        
        % ======================================================================
        %> @brief Build the total scattering source for a group.
        %
        %> In some cases, including all scattering is required, as is the case
        %> when performing multigroup Krylov solves.  This constructs
        %> @f[
        %>     q_g = \sum^G_{g'} \mathbf{S}_{gg'}\phi_{g'} \, .
        %> @f]
        %>  
        %> @param   g       Group for this problem.  (I.e. row in MG).
        %> @param   phi     Multigroup fluxes.
        %> @return          Total scatter source.        
        % ======================================================================
        function q = build_total_scatter_source(this, g, phi)
            q = 0.0;
            for gp = lower(this.d_mat, g):upper(this.d_mat, g)
                q = q + phi(:, gp) .* this.d_sigma_s{g}(:, gp);
            end        
        end % end function build_all_scatter_source   
        
        
    end
    
end