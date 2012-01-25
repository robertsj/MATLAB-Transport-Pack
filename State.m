%> @file  State.m
%> @brief State class definition.
% ==============================================================================
%> @brief Represents the state of the problem.
%
%> This contains everything that needs to be known throughout most of the
%> problem.  The scalar flux (and maybe moments later on), an eigenvalue, and
%> boundary fluxes all live here.
% ==============================================================================
classdef State < handle

   properties (Access = private)
       %> Cell center scalar flux (only isotropic at this point)
       d_phi
       %> Eigenvalue
       d_eigenvalue = 1.0;
       %> Mesh
       d_mesh
       %> Boundary fluxes.
       d_boundary
       %> Number of groups.
       d_number_groups
   end
   
   methods
       
        % ======================================================================
        %> @brief Class constructor
        %>
        %> More detailed description of what the constructor does.
        %>
        %> @param input         User input.
        %> @param mesh          Problem mesh.
        %> @param quadrature    Angular quadrature.
        %>
        %> @return Instance of the State class.
        % ======================================================================
        function this = State(input, mesh, quadrature)
            this.d_number_groups = get(input, 'number_groups');
            this.d_mesh = mesh;
            % Here, the number of cells is either the fine mesh count for
            % discrete ordinates or the flat source regions in MOC.
            if meshed(mesh)
                n = number_cells(mesh);
            else
                n = number_regions(mesh);
            end
            this.d_phi = zeros(n, this.d_number_groups);
            
        end
        
        % ======================================================================
        %> @brief Get a group flux vector.
        %
        %> @param  g    Energy group.
        %> @return      Flux vector for all cells.
        % ======================================================================
        function f = flux(this, g)
            f = this.d_phi(:, g);
        end
        
        % ======================================================================
        %> @brief Get the last computed eigenvalue.
        %
        %> @return  Eigenvalue (i.e. keff)
        % ======================================================================
        function k = eigenvalue(this)
            k = this.d_eigenvalue;
        end
        
        % ======================================================================
        %> @brief Get the boundary fluxes.
        %
        %> @return  Eigenvalue (i.e. keff)
        % ======================================================================
        function b = boundary(this)
            b = this.d_boundary; 
        end
        
        % ======================================================================
        %> @brief Set a group flux vector.
        %
        %> @param  phi  Flux vector.
        %> @param  g    Energy group.
        % ======================================================================
        function this = set_phi(this, phi, g)
            this.d_phi(:, g) = phi(:);
        end
        
        % ======================================================================
        %> @brief Set the eigenvalue.
        %
        %> @param  eigenvalue  Eigenvalue estimate.
        % ======================================================================
        function this = set_eigenvalue(this, eigenvalue)
            this.d_eigenvalue = eigenvalue;
        end
        
        % ======================================================================
        %> @brief Set the boundary fluxes.
        %
        %> @param  boundary  Boundary flux object.
        % ======================================================================
        function this = set_boundary(this, boundary)
            this.d_boundary = boundary;
        end        
        
   end
    
end