%> @file  Equation.m
%> @brief Equation class definition.
% ==============================================================================
%> @brief Base equation class for spatial discretization.
%
%> More...
% ==============================================================================
classdef Equation < handle
    
    properties (Constant)
       HORZ = 1;
       VERT = 2;
    end
    
    properties (Access = protected)   
        %> Problem mesh
        d_mesh
        %> Material definitions
        d_mat
        %> Quadrature
        d_quadrature
        %> Current mu value
        d_mu
        %> current eta value
        d_eta
        %> Current ksi value
        d_ksi
        %> Weighted diamond difference parameter. 
        d_alpha
        %> Material map
        d_mat_map
    end
    
    methods
        
        % ======================================================================
        %> @brief Class constructor
        %>
        %> Set the mesh and material.
        %>
        %> @param mesh        	Problem mesh.
        %> @param mat          	Material definitions.
        %>
        %> @return Instance of the Equation class.
        % ======================================================================
        function obj = Equation(mesh, mat, quadrature)
            obj.d_mesh = mesh;
            obj.d_mat  = mat;
            obj.d_quadrature = quadrature;
            if meshed(mesh)
                obj.d_mat_map = mesh_map(mesh, 'MATERIAL');
            end
        end
        
        % ======================================================================
        %> @brief Setup the equations for a group.
        %>
        %> @param group     Current group.
        % ======================================================================
        obj = setup_group(obj, group)     
        
        % ======================================================================
        %> @brief Setup the equations for an octant.
        %>
        %> @param octant    Current octant.
        % ======================================================================
        obj = setup_octant(obj, octant)
        
        % ======================================================================
        %> @brief Setup the equations for an angle.
        %>
        %> @param mu    Cosine with respect to x axis.
        %> @param eta   Cosine with respect to y axis.
        %> @param xi    Cosine with respect to z axis.
        % ======================================================================
        obj = setup_angle(obj, mu, eta, xi)
        
        % ======================================================================
        %> @brief Solve for the cell-center and outgoing edge fluxes.
        %
        %> @param psi_in    Incident flux vector, [horz, vert]
        %> @param s         Cell source
        %> @param i         Cell x index
        %> @param j         Cell y index
        %> @param g         Group index
        %>
        %> @return Cell center angular flux and outgoing edge fluxes.
        % ======================================================================
        [psi_out psi_center] = solve(obj, g, psi_in, s, i, j, k)
        
        % ======================================================================
        %> @brief Number of unknowns per cell.
        %
        %> This helps to allow various equation types to be used
        %> generically, without knowing the number of unknowns per cell.
        %> For the simplest schems (e.g. diamond difference), there is only
        %> one unknown per cell.  For other schemes (e.g. linear
        %> discontinuous, there will be more.
        %>
        %> @return Number of unknowns.
        % ======================================================================
        u = number_unknowns(obj);
           
    end
    
end