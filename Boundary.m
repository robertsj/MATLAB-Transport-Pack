%> @file  Boundary.m
%> @brief Boundary class definition.
% ==============================================================================
%> @brief Holds the incident and outgoing boundary fluxes.
%
%> More here...
% ==============================================================================
classdef Boundary < handle
    
    properties (Constant)
       IN   =  1;
       OUT  = -1;
    end

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
            
            
            ng = get(input, 'number_groups');
            % Left and right
            obj.d_boundary_flux{mesh.LEFT} = zeros(dim(mesh, 2), ...
                number_angles(quadrature), ng);
            obj.d_boundary_flux{mesh.RIGHT} = zeros(dim(mesh, 2), ...
                number_angles(quadrature), ng);     
            
            % Top and bottom
            obj.d_boundary_flux{mesh.TOP} = zeros(dim(mesh, 1), ...
                number_angles(quadrature), ng);
            obj.d_boundary_flux{mesh.BOTTOM} = zeros(dim(mesh, 1), ...
                number_angles(quadrature), ng);  
           
            % Boundary conditions
            obj.d_bc = cell(2*mesh.DIM);
            
            if strcmp(get(input, 'bc_left'), 'vacuum')
                obj.d_bc{mesh.LEFT}   = ...
                    Vacuum(obj, mesh, quadrature, mesh.LEFT);
            elseif strcmp(get(input, 'bc_left'), 'reflect')
                obj.d_bc{mesh.LEFT}   = ...
                    Reflective(obj, mesh, quadrature, mesh.LEFT);
            end
            if strcmp(get(input, 'bc_right'), 'vacuum')
                obj.d_bc{mesh.RIGHT}   = ...
                    Vacuum(obj, mesh, quadrature, mesh.RIGHT);
            elseif strcmp(get(input, 'bc_right'), 'reflect')
                obj.d_bc{mesh.RIGHT}   = ...
                    Reflective(obj, mesh, quadrature, mesh.RIGHT);
            end   
            if mesh.DIM > 1
                if strcmp(get(input, 'bc_bottom'), 'vacuum')
                    obj.d_bc{mesh.BOTTOM}   = ...
                        Vacuum(obj, mesh, quadrature, mesh.BOTTOM);
                elseif strcmp(get(input, 'bc_bottom'), 'reflect')
                    obj.d_bc{mesh.BOTTOM}   = ...
                        Reflective(obj, mesh, quadrature, mesh.BOTTOM);
                end
                if strcmp(get(input, 'bc_top'), 'vacuum')
                    obj.d_bc{mesh.TOP}   = ...
                        Vacuum(obj, mesh, quadrature, mesh.TOP);
                elseif strcmp(get(input, 'bc_top'), 'reflect')
                    obj.d_bc{mesh.TOP}   = ...
                        Reflective(obj, mesh, quadrature, mesh.TOP);
                end
            end
            
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
            for side = 1:2*obj.d_mesh.DIM;
                update(obj.d_bc{side});
            end
        end
        
        % Getters
        
        
        % ======================================================================
        %> @brief Get a vertical face incident/exiting boundary flux.
        %
        %> @param o     Octant index.
        %> @param a     Angle index (within an octant).
        %> @param inout Incoming/outgoing switch 
        %> @return      Boundary flux array.
        % ======================================================================
        % example:  We are starting from the left, octant 1, going right.  We
        % call this function with o=1, inout=IN=1.  We get LEFT.
        %
        function f = get_psi_v(obj, o, a, inout)
            if obj.d_quadrature.octant(o, 1) == inout % 1/-1
                side = obj.d_mesh.LEFT;     % We are starting from the left.
            else
                side = obj.d_mesh.RIGHT;    % We are starting from the right.
            end
            %disp(['getting in(out)',num2str(inout),' for side=', num2str(side)]);
            % Cardinal angular index.
            angle = index(obj.d_quadrature, o, a);
            f = obj.d_boundary_flux{side}(:, angle, obj.d_g);           
        end
        
        % ======================================================================
        %> @brief Get a horizontal face incident/exiting boundary flux.
        %
        %> @param o     Octant index.
        %> @param a     Angle index (within an octant).
        %> @param inout Incoming/outgoing switch         
        %> @return      Boundary flux array.
        % ======================================================================
        function f = get_psi_h(obj, o, a, inout)
            if obj.d_quadrature.octant(o, 2) == inout % 1/-1
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
        %> @param inout Incoming/outgoing switch         
        % ======================================================================
        % example: sweeping left to right, i.e octant 1.  At the end, we set the
        % right boundary octant 1 flux.  inout=OUT=-1.  octant(1) = 1.  Hence,
        % the logic below finds RIGHT.  
        %
        function obj = set_psi_v(obj, o, a, f, inout)
            if obj.d_quadrature.octant(o, 1) == inout  
                side = obj.d_mesh.LEFT;     % We are exiting the left.
            else
                side = obj.d_mesh.RIGHT;    % We are exiting the right.
            end
           % disp(['setting in(out)',num2str(inout),' for side=', num2str(side)]);
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
        %> @param inout Incoming/outgoing switch         
        % ======================================================================
        function obj = set_psi_h(obj, o, a, f, inout)
            if obj.d_quadrature.octant(o, 2) == inout
                side = obj.d_mesh.BOTTOM;   % Entering/exiting the bottom.
            else
                side = obj.d_mesh.TOP;      % Entering/exiting the top.
            end
            % Cardinal angular index.
            angle = index(obj.d_quadrature, o, a);
            obj.d_boundary_flux{side}(:, angle, obj.d_g) = f;      
        end   
        
    end
 
end