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
            d  = get(input, 'dimension');
            na = number_angles(quadrature);
            this.d_boundary_flux = cell(2*d, 1);
                        
            % Initialize boundary fluxes.
            if (d == 1)
                this.d_boundary_flux{mesh.LEFT}   = zeros(1, na, ng);
                this.d_boundary_flux{mesh.RIGHT}  = zeros(1, na, ng);
            elseif (d == 2)
                this.d_boundary_flux{mesh.LEFT}   = zeros(dim(mesh, 2), na, ng);
                this.d_boundary_flux{mesh.RIGHT}  = zeros(dim(mesh, 2), na, ng);
                this.d_boundary_flux{mesh.TOP}    = zeros(dim(mesh, 1), na, ng);
                this.d_boundary_flux{mesh.BOTTOM} = zeros(dim(mesh, 1), na, ng);
            elseif (d == 3)
                this.d_boundary_flux{mesh.LEFT}   = ...
                    zeros(dim(mesh, 2)*dim(mesh, 3), na, ng);
                this.d_boundary_flux{mesh.RIGHT}  = ...
                    zeros(dim(mesh, 2)*dim(mesh, 3), na, ng);
                this.d_boundary_flux{mesh.TOP}    = ...
                    zeros(dim(mesh, 1)*dim(mesh, 3), na, ng);
                this.d_boundary_flux{mesh.BOTTOM} = ... 
                    zeros(dim(mesh, 2)*dim(mesh, 3), na, ng);
            else
                error('user:input', 'Invalid dimension');
            end
            
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
            for side = 1:2*d;
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
        
        % ======================================================================
        %> @brief Get total leakage in current group.
        % ======================================================================
        function leak = get_leakage(this, inout)
            if nargin == 1
                inout = Boundary.OUT;
            end
            w   = weight_octant(this.d_quadrature);
            leak = zeros(1, 4); 
            for o = 1:4
                if this.d_quadrature.octant(o, 1) == inout % 1/-1
                    sidev = this.d_mesh.LEFT;     % We are starting from the left.
                else
                    sidev = this.d_mesh.RIGHT;    % We are starting from the right.
                end                
                if this.d_quadrature.octant(o, 2) == inout % 1/-1
                    sideh = this.d_mesh.BOTTOM;   % We are starting from the bottom.
                else
                    sideh = this.d_mesh.TOP;      % We are starting from the top.
                end                
                width = widths(this.d_mesh); 
                wx = width{1};
                wy = width{2};
                psi_v = get_psi_v_octant(this, o, inout);
                psi_h = get_psi_h_octant(this, o, inout);
                mu  = angle_octant(this.d_quadrature, 1);
                eta = angle_octant(this.d_quadrature, 2);
                w   = weight_octant(this.d_quadrature);
                % leak = dx(n,1) .* psi(n,m)*muw(m, 1)
                leak(sidev) = leak(sidev) + wy' * psi_v*(mu.*w);
                leak(sideh) = leak(sideh) + wx' * psi_h*(eta.*w);
                %
%                 for a = 1:number_angles_octant(this.d_quadrature)
%                     [mu, eta] = angle(this.d_quadrature, o, a);
%                     psi_v = get_psi_v(this, o, a, inout);
%                     for j = 1:length(psi_v(:, 1))
%                         leak(sidev) = leak(sidev) + dy(this.d_mesh, j) * psi_v(j)*abs(mu)*w(a);
%                     end
%                     psi_h = get_psi_h(this, o, a, inout);
%                     for i = 1:length(psi_h(:, 1))
%                         leak(sideh) = leak(sideh) + dx(this.d_mesh, i) * psi_h(i)*abs(eta)*w(a);
%                     end
%                 end
            end
            
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
            try
            this.d_boundary_flux{side}(:, ...
                this.d_octant_indices(o, 1):this.d_octant_indices(o, 2), ...
                this.d_g) = f(:, :);   
            catch ME
                disp('lala')
            end
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