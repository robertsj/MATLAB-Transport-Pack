%> @file  SweepBase.m
%> @brief SweepBase class definition.
% ==============================================================================
%> @brief Base sweep class for discrete ordinates.
%
%> The within-group transport equation is 
%> \f[
%>      \mathbf{T}\psi = Q \, ,
%> \f]
%> where \f$ \mathbf{T} \f$ is the streaming and collision operator and 
%> \f$ Q \f$ is a discrete representation of all source contributions.
%>
%> To invert the operator \f$ \mathbf{T} \f$, we "sweep" over the mesh for all 
%> angles,
%> which gives us updated angular fluxes in each cell.  However, we don't store
%> the angular flux, but rather add its contribution to the scalar flux 
%> directly since storing the angular flux is too expensive.
% ==============================================================================
classdef SweepBase < handle

    properties (Access = public)
        
        %> User input.
        d_input      
        %> Cartesian mesh.
        d_mesh                  
        %> Materials.
        d_mat                   
        %> Angular mesh.
        d_quadrature            
        %> Spatial discretization
        d_equation                      
        %> Boundary fluxes
        d_boundary
        
        %> Number of cells in the x direction
        d_nx
        %> Number of cells in the y direction
        d_ny
        %> Number of cells in the z direction
        d_nz
        
        %> Number of sweeps
        d_number_sweeps;
        
    end
    
    methods
        
        % ======================================================================
        %> @brief Class constructor
        %
        %> @param   input       Input database.
        %> @param   mesh      	Cartesian mesh.
        %> @param   mat       	Material definitions.
        %> @param   quadrature  Angular mesh.
        %> @param   boundary    Boundary flux container.
        %> @param   equation    Discretization
        %>
        %> @return              Instance of the SweepBase class.
        % ======================================================================
        function obj = SweepBase(input, mesh, mat, quadrature, boundary, ...
                                 equation)
            
            obj.d_input      = input;
            obj.d_mesh       = mesh;
            obj.d_mat        = mat;
            obj.d_quadrature = quadrature;
            obj.d_boundary   = boundary;
            obj.d_equation   = equation;
            
            % Store some things from the mesh.
            obj.d_nx = number_cells_x(mesh);
            obj.d_ny = number_cells_y(mesh);
            obj.d_nz = number_cells_z(mesh);
            
            obj.d_number_sweeps = 0;

        end
        
        % ======================================================================
        %> @brief Sweep the mesh for all angles.
        %
        %> This performs the action of \f$ \mathbf{T}^{-1} \f$ on a given
        %> discrete right hand side vector.  Currently, a single vector
        %> applicable to all angles is provided, since we work only with 
        %> isotropic sources and scattering.
        %>
        %> @param s         Discrete sweep source.
        %> @param g         Group of this problem.
        %> @return          Updated group flux.
        % ======================================================================
        phi = sweep(obj, s, g)

        %> @name Getters
        %> @{
        
        function n = number_sweeps(this)
            n = this.d_number_sweeps;
        end
        
        %> @}
        
        
    end
   
    
    
    
end