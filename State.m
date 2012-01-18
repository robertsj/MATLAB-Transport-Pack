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
       %> Outgoing boundary fluxes.  This is an array of cells, one per surface.
       d_boundary_flux
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
        function obj = State(input, mesh, quadrature)
            obj.d_mesh = mesh;
            % MATLAB uses column-major ordering, so in many cases, it might be
            % best to order space as fastest changing index
            obj.d_phi = zeros(number_cells(mesh), input.number_groups);              
        end
        
        function f = flux(obj, g)
            % Returns a group flux for all cells.
            f = obj.d_phi(:, g);
        end
        
        function k = eigenvalue(obj)
            % Returns the current eigenvalue.
            k = obj.d_eigenvalue;
        end
        
        function obj = set_eigenvalue(obj, eigenvalue)
            obj.d_eigenvalue = eigenvalue;
        end
        
        function obj = set_phi(obj, phi, g)
            obj.d_phi(:, g) = phi(:);
        end
        
        function obj = set_boundary_flux(obj, dim, psi, a, g)
            obj.d_boundary_flux{dim}(:, a, g) = psi;
        end
        
   end
    
end