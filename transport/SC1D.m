%> @file  SC1D.m
%> @brief SC1D class definition.
% ==============================================================================
%> @brief Step characteristic approximation in one dimension.
%
%> 
%> Finish me.
% ==============================================================================
classdef SC1D < Equation
    
    properties
        d_con_x
        d_sig
    end
    
    properties (Constant)
%         HORZ = 1;
%         VERT = 2;
    end
    
    methods
        
        % ======================================================================
        %> @brief Class constructor
        %
        %> Set the mesh and material.
        %>
        %> @param mesh        	Problem mesh.
        %> @param mat          	Material definitions.
        %>
        %> @return Instance of the SC1D class.
        % ======================================================================
        function obj = SC1D(mesh, mat)
            % Call the base class.
            obj = obj@Equation(mesh, mat); 
            % Presize coefficient vectors.
            obj.d_con_x = zeros(number_cells_x(mesh), 1);
            obj.d_sig   = zeros(number_cells_x(mesh), 1);
        end
        
        % ======================================================================
        %> @brief Setup the equations for a group.
        %
        %> Here, we'll go through the grid and produce a fine mesh matrix
        %> of total cross-sections.  This isn't the best thing for memory,
        %> but it cuts down a lot on the time within solve.
        %>
        %> @param group     Current group.
        % ======================================================================
        function obj = setup_group(obj, group)
            for i = 1:number_cells_x(obj.d_mesh)
                obj.d_sig(i) = sigma_t(obj.d_mat, obj.d_mat_map(i), group);
            end
        end
        
        % ======================================================================
        %> @brief Setup the equations for an octant.
        %>
        %> @param octant    Current octant.
        % ======================================================================
        function obj = setup_octant(obj, octant)
            % Nothing here for now.   
        end
        
        % ======================================================================
        %> @brief Setup the equations for an angle.
        %>
        %> @param mu    Cosine with respect to x axis.
        %> @param eta   Cosine with respect to y axis.
        % ======================================================================
        function obj = setup_angle(obj, mu, eta)
            % Get the widths from mesh.
            w = widths(obj.d_mesh);
            % Build the two constants.
            obj.d_con_x = 2*abs(mu)./w{1}; 
        end
        
        % ======================================================================
        %> @brief Solve for the cell-center and outgoing edge fluxes.
        %>
        %> @param g         Group index
        %> @param psi_in    Incident flux vector
        %> @param s         Cell source
        %> @param i         Cell x index
        %> @return Cell center angular flux and outgoing edge fluxes.
        % ======================================================================
        function [psi_out, psi_center] = solve(obj, g, psi_in, s, i)
            coef = 1.0 / (obj.d_sig(i) + obj.d_con_x(i));
            psi_center = coef * (s + obj.d_con_x(i) * psi_in(1));
            psi_out(1) = 2*psi_center - psi_in(1);
        end
        
    end
    
end