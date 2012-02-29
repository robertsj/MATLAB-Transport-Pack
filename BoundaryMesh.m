%> @file  BoundaryMesh.m
%> @brief BoundaryMesh class definition.
% ==============================================================================
%> @brief Holds the incident and outgoing boundary fluxes for a mesh.
%>
%> For now, this class is largely limited to 1-d and 2-d problems.  Access is
%> with respect to "vertical" (1-d/2-d) and "horizontal (2-d only) edges,
%> specified as a function of octant and whether the incident or outgoing flux 
%> is desired.
% ==============================================================================
classdef BoundaryMesh < Boundary
    
    properties
        %> Octant indices
        d_octant_indices
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
        %> @return Instance of the BoundaryMesh class.
        % ======================================================================
        function this = BoundaryMesh(input, mesh, quadrature)
        
            % Call the base class.
            this = this@Boundary(input, mesh, quadrature);
            
            % Initialize boundary flux vectors.
            ng = get(input, 'number_groups');
            this.d_boundary_flux = cell(2*mesh.DIM, 1);
            % Left and right
            this.d_boundary_flux{mesh.LEFT} = zeros(dim(mesh, 2), ...
                number_angles(quadrature), ng);
            this.d_boundary_flux{mesh.RIGHT} = zeros(dim(mesh, 2), ...
                number_angles(quadrature), ng);     
            
            % Top and bottom
            this.d_boundary_flux{mesh.TOP} = zeros(dim(mesh, 1), ...
                number_angles(quadrature), ng);
            this.d_boundary_flux{mesh.BOTTOM} = zeros(dim(mesh, 1), ...
                number_angles(quadrature), ng);  
            
            % Build octant indices
            this.d_octant_indices = ...
                zeros(this.d_quadrature.number_octants(), 2);
            
            for i = 1:this.d_quadrature.number_octants()
                 this.d_octant_indices(i, 1) = ...
                     1 + (i-1)*this.d_quadrature.number_angles_octant();
                 this.d_octant_indices(i, 2) = ...
                     i*this.d_quadrature.number_angles_octant();
            end
            
            % Initialize the boundary conditions
            for side = 1:2*this.d_mesh.DIM;
                initialize(this.d_bc{side});
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
        function f = get_psi_v(this, o, a, inout)
            if this.d_quadrature.octant(o, 1) == inout % 1/-1
                side = this.d_mesh.LEFT;     % We are starting from the left.
            else
                side = this.d_mesh.RIGHT;    % We are starting from the right.
            end
            %disp(['getting in(out)',num2str(inout),' for side=', num2str(side)]);
            % Cardinal angular index.
            angle = index(this.d_quadrature, o, a);
            f = this.d_boundary_flux{side}(:, angle, this.d_g);           
        end
        
        % ======================================================================
        %> @brief Get a horizontal face incident/exiting boundary flux.
        %
        %> @param o     Octant index.
        %> @param a     Angle index (within an octant).
        %> @param inout Incoming/outgoing switch         
        %> @return      Boundary flux array.
        % ======================================================================
        function f = get_psi_h(this, o, a, inout)
            if this.d_quadrature.octant(o, 2) == inout % 1/-1
                side = this.d_mesh.BOTTOM;   % We are starting from the bottom.
            else
                side = this.d_mesh.TOP;      % We are starting from the top.
            end
            % Cardinal angular index.
            angle = index(this.d_quadrature, o, a);
            f = this.d_boundary_flux{side}(:, angle, this.d_g);         
        end
        
        % ======================================================================
        %> @brief Get a vertical face incident/exiting boundary flux.
        %
        %> @param o     Octant index.
        %> @param a     Angle index (within an octant).
        %> @param inout Incoming/outgoing switch 
        %> @return      Boundary flux array.
        % ======================================================================
        function f = get_psi_v_octant(this, o, inout)
            if this.d_quadrature.octant(o, 1) == inout % 1/-1
                side = this.d_mesh.LEFT;     % We are starting from the left.
            else
                side = this.d_mesh.RIGHT;    % We are starting from the right.
            end
            f = this.d_boundary_flux{side}(:, ...
                this.d_octant_indices(o, 1):this.d_octant_indices(o, 2), ...
                this.d_g);           
        end
        
        % ======================================================================
        %> @brief Get a horizontal face incident/exiting boundary flux.
        %
        %> @param o     Octant index.
        %> @param a     Angle index (within an octant).
        %> @param inout Incoming/outgoing switch         
        %> @return      Boundary flux array.
        % ======================================================================
        function f = get_psi_h_octant(this, o, inout)
            if this.d_quadrature.octant(o, 2) == inout % 1/-1
                side = this.d_mesh.BOTTOM;   % We are starting from the bottom.
            else
                side = this.d_mesh.TOP;      % We are starting from the top.
            end
            % Cardinal angular index.      
            f = this.d_boundary_flux{side}(:, ...
                this.d_octant_indices(o, 1):this.d_octant_indices(o, 2), ...
                this.d_g);              
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
        function this = set_psi_v(this, o, a, f, inout)
            if this.d_quadrature.octant(o, 1) == inout  
                side = this.d_mesh.LEFT;     % We are exiting the left.
            else
                side = this.d_mesh.RIGHT;    % We are exiting the right.
            end
            % Cardinal angular index.
            angle = index(this.d_quadrature, o, a);
            this.d_boundary_flux{side}(:, angle, this.d_g) = f;
        end
        
        % ======================================================================
        %> @brief Set a horizontal face incident boundary flux.
        %
        %> @param o     Octant index.
        %> @param a     Angle index (within an octant).
        %> @param f     Flux array.
        %> @param inout Incoming/outgoing switch         
        % ======================================================================
        function this = set_psi_h(this, o, a, f, inout)
            if this.d_quadrature.octant(o, 2) == inout
                side = this.d_mesh.BOTTOM;   % Entering/exiting the bottom.
            else
                side = this.d_mesh.TOP;      % Entering/exiting the top.
            end
            % Cardinal angular index.
            angle = index(this.d_quadrature, o, a);
            this.d_boundary_flux{side}(:, angle, this.d_g) = f;      
        end   
        
        % ======================================================================
        %> @brief Set a vertical face incident boundary flux for all angles.      
        %
        %> @param o     Octant index.
        %> @param f     Flux array.
        %> @param inout Incoming/outgoing switch         
        % ======================================================================
        function this = set_psi_v_octant(this, o, f, inout)
            if this.d_quadrature.octant(o, 1) == inout  
                side = this.d_mesh.LEFT;     % We are exiting the left.
            else
                side = this.d_mesh.RIGHT;    % We are exiting the right.
            end
            this.d_boundary_flux{side}(:, ...
                this.d_octant_indices(o, 1):this.d_octant_indices(o, 2), ...
                this.d_g) = f(:, :);      
        end
        
        % ======================================================================
        %> @brief Set a horizontal face incident boundary flux for all angles. 
        %
        %> @param o     Octant index.
        %> @param f     Flux array.
        %> @param inout Incoming/outgoing switch         
        % ======================================================================
        function this = set_psi_h_octant(this, o, f, inout)
            if this.d_quadrature.octant(o, 2) == inout
                side = this.d_mesh.BOTTOM;   % Entering/exiting the bottom.
            else
                side = this.d_mesh.TOP;      % Entering/exiting the top.
            end
            this.d_boundary_flux{side}(:, ...
                this.d_octant_indices(o, 1):this.d_octant_indices(o, 2), ...
                this.d_g) = f(:, :);      
        end    
    end
    
end