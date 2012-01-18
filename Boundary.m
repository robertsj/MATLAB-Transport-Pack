%> @file  Boundary.m
%> @brief Boundary class definition.
% ==============================================================================
%> @brief Holds the incident and outgoing boundary fluxes.
%
%> More here...
% ==============================================================================
classdef Boundary < handle

    properties (Access = private)
        %> Spatial mesh.
        d_mesh
        %> Angular mesh.
        d_quadrature
        %> Boundary flux cell array
        d_boundary_flux
        %> Cell array of \ref BoundaryCondition objects for each surface
        d_bc  
        %> Current group.
        d_g = 0;
    end
    
    methods
        
        % ======================================================================
        %> @brief Class constructor
        %>
        %> All boundary flux arrays are sized, and the boundary conditions are
        %> constructed for each surface.
        %>
        %> @param input         User input.
        %> @param mesh          Spatial mesh.
        %> @param quadrature  	Angular mesh.
        %>
        %> @return Instance of the Boundary class.
        % ======================================================================
        function obj = Boundary(input, mesh, quadrature)
        
            obj.d_mesh = mesh;
            obj.d_quadrature = quadrature;
            obj.d_boundary_flux = cell(2*mesh.DIM, 1);
            
            % Left and right
            obj.d_boundary_flux{mesh.LEFT} = zeros(dim(mesh, 2), ...
                number_angles(quadrature), input.number_groups);
            obj.d_boundary_flux{mesh.RIGHT} = zeros(dim(mesh, 2), ...
                number_angles(quadrature), input.number_groups);     
            
            % Top and bottom
            obj.d_boundary_flux{mesh.TOP} = zeros(dim(mesh, 1), ...
                number_angles(quadrature), input.number_groups);
            obj.d_boundary_flux{mesh.BOTTOM} = zeros(dim(mesh, 1), ...
                number_angles(quadrature), input.number_groups);  
           
            % Boundary conditions
            obj.d_bc = cell(2*mesh.DIM);
            obj.d_bc{mesh.LEFT}   = Vacuum(obj, mesh.LEFT);
            obj.d_bc{mesh.RIGHT}  = Vacuum(obj, mesh.RIGHT);
            obj.d_bc{mesh.BOTTOM} = Vacuum(obj, mesh.BOTTOM);
            obj.d_bc{mesh.TOP}    = Vacuum(obj, mesh.TOP);
            
        end
        
        % ======================================================================
        %> @brief Initialize for a within-group solve.
        %>
        %> More detailed description of what the constructor does.
        %>
        %> @param input         User input.
        %> @param mesh          Spatial mesh.
        %> @param quadrature  	Angular mesh.
        %>
        %> @return Instance of the Boundary class.
        % ======================================================================
        function obj = initialize(obj, g)
            obj.d_g = g;
        end
        
        
        % ======================================================================
        %> @brief Set the incident boundary fluxes.
        %>
        %> This is called after every sweep.  The sweeper updates all the
        %> outgoing angular fluxes.  This routine updates the incident fluxes
        %> according to the boundary condition.
        %>
        %> @param order         Quadrature order. This differs from quadrature 
        %>                      to quadrature, but e.g. for level symmetric, 
        %>                      it's the number of unique directional cosines.
        %>
        %> @return Instance of the Quadrature class.
        % ======================================================================  
        function obj = set(obj)
            for side = 1:4
                update(obj.d_bc{side});
            end
        end
        
        % Getters
        
        
        % ======================================================================
        %> @brief Get a vertical face incident boundary flux.
        %
        %> @param o     Octant index.
        %> @param a     Angle index (within an octant).
        %> @return      Boundary flux array.
        % ======================================================================
        function f = get_psi_v(obj, o, a)
            if obj.d_quadrature.octant(o, 1) == 1
                side = obj.d_mesh.LEFT;     % We are starting from the left.
            else
                side = obj.d_mesh.RIGHT;    % We are starting from the right.
            end
            % Cardinal angular index.
            angle = index(obj.d_quadrature, o, a);
            f = obj.d_boundary_flux{side}(:, angle, obj.d_g);         
        end
        
        % ======================================================================
        %> @brief Get a horizontal face incident boundary flux.
        %
        %> @param o     Octant index.
        %> @param a     Angle index (within an octant).
        %> @return      Boundary flux array.
        % ======================================================================
        function f = get_psi_h(obj, o, a)
            if obj.d_quadrature.octant(o, 2) == 1
                side = obj.d_mesh.BOTTOM;   % We are starting from the bottom.
            else
                side = obj.d_mesh.TOP;      % We are starting from the top.
            end
            % Cardinal angular index.
            angle = index(obj.d_quadrature, o, a);
            f = obj.d_boundary_flux{side}(:, angle, obj.d_g);         
        end
        
        % Setters
        
        % ======================================================================
        %> @brief Set a vertical face incident boundary flux.
        %
        %> @param o     Octant index.
        %> @param a     Angle index (within an octant).
        %> @param f     Flux array.
        % ======================================================================
        function obj = set_psi_v(obj, o, a, f)
            if obj.d_quadrature.octant(o, 1) == -1
                side = obj.d_mesh.LEFT;     % We are exiting the left.
            else
                side = obj.d_mesh.RIGHT;    % We are exiting the right.
            end
            % Cardinal angular index.
            angle = index(obj.d_quadrature, o, a);
            obj.d_boundary_flux{side}(:, angle, obj.d_g) = f;
        end
        
        % ======================================================================
        %> @brief Set a horizontal face incident boundary flux.
        %
        %> @param o     Octant index.
        %> @param a     Angle index (within an octant).
        %> @param f     Flux array.
        % ======================================================================
        function obj = set_psi_h(obj, o, a, f)
            if obj.d_quadrature.octant(o, 2) == -1
                side = obj.d_mesh.BOTTOM;   % We are exiting the bottom.
            else
                side = obj.d_mesh.TOP;      % We are exiting the top.
            end
            % Cardinal angular index.
            angle = index(obj.d_quadrature, o, a);
            obj.d_boundary_flux{side}(:, angle, obj.d_g) = f;      
        end   
        
    end
 
end