%> @file  BoundaryCondition.m
%> @brief Boundary condition class definition.
% ==============================================================================
%> @brief Base boundary condition class.
% ==============================================================================
classdef BoundaryCondition < handle

    properties (Access = protected) 
        %>
        d_input
        %> Boundary flux thisect.
        d_boundary
        %> Mesh 
        d_mesh
        %> Quadrature
        d_quadrature
        % Surface identifier.
        d_side = 0 
    end
    
    methods
        
        % ======================================================================
        %> @brief Class constructor
        %>
        %> More detailed description of what the constructor does.
        %>
        %> @param boundary      Boundary flux class.
        %> @param side          Surface identifier.
        %>
        %> @return Instance of the Boundary class.
        % ======================================================================
        function this = BoundaryCondition(boundary, input, mesh, ...
                                          quadrature, side)
            this.d_input         = input;
            this.d_boundary      = boundary;
            this.d_mesh          = mesh;
            this.d_quadrature    = quadrature;
            this.d_side          = side;
        end

        % ======================================================================
        %> @brief Initialize the boundary condition.
        %>
        %> For fixed boundaries, this is the only relevant call. Other
        %> need not implement this.
        % ======================================================================        
        function this = initialize(this)
            
        end        
        
        % ======================================================================
        %> @brief Update the boundary flux.
        % ======================================================================
        this = update(this);
        
    end
 
end