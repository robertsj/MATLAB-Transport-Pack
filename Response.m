%> @file  Response.m
%> @brief Response boundary condition class definition.
% ==============================================================================
%> @brief Reflective boundary condition.
%
%> 
% ==============================================================================
classdef Response < BoundaryCondition
    
    properties (Access = protected)
        %> Incident octants
        d_octants
        %> Am I set? (Default false)
        d_set = 0
        %> Group
        d_order_group
        %> Spatial order
        d_order_space
        %> Polar order
        d_order_polar 
        %> Azimuthal order
        d_order_azimuth
        %> Spatial basis
        d_basis_space
        %> Polar angle basis        
        d_basis_polar
        %> Azimuthal angle basis        
        d_basis_azimuth
        %> Number of spatial cells   
        d_number_space
        %> Number of polar angles   
        d_number_polar
        %> Number of azimuthal angles
        d_number_azimuth
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
        %> @return Instance of the Response class.
        % ======================================================================
        function this = Response(boundary, input, mesh, quadrature, side)
            this = this@BoundaryCondition(boundary, input, mesh, quadrature, side);

            % Build the bases, normalized DLP's.   
            this.d_order_group   = get(input, 'rf_order_group');
            this.d_order_space   = get(input, 'rf_order_space');
            this.d_order_polar   = get(input, 'rf_order_polar');            
            this.d_order_azimuth = get(input, 'rf_order_azimuth');
            
            
            % ======================
            % 1-D
            % ======================
            if mesh.DIM == 1
                if side == Mesh.LEFT
                    this.d_octants(1, 1) =  1; % I go in octant one, the right
                end
                if side == Mesh.RIGHT
                    this.d_octants(1, 1) =  2; % I go in octant two, the left
                end
                
                % Building all orders for now.
                n = number_angles_octant(quadrature);
                this.d_basis_polar = DiscreteLP(n-1);
           
            % ======================
            % 2-D
            % ======================         
            elseif mesh.DIM == 2
                if side == Mesh.LEFT
                    this.d_octants(1, 1) =  1; % my "left"
                    this.d_octants(2, 1) =  4; % my "right"
                elseif side == Mesh.RIGHT
                    this.d_octants(1, 1) =  3; 
                    this.d_octants(2, 1) =  2;       
                elseif side == Mesh.BOTTOM
                    this.d_octants(1, 1) =  2;
                    this.d_octants(2, 1) =  1; 
                elseif side == Mesh.TOP
                    this.d_octants(1, 1) =  4; 
                    this.d_octants(2, 1) =  3; 
                else
                    error('Wrong side for 2-D Reflective')
                end
                
                % Building all orders for now.
                if side < Mesh.BOTTOM
                    this.d_number_space = number_cells_x(mesh);
                else
                    this.d_number_space = number_cells_y(mesh);
                end
                this.d_number_polar   = get(input, 'quad_number_polar');
                this.d_number_azimuth = get(input, 'quad_number_azimuth');                
                
                this.d_basis_space   = DiscreteLP(this.d_number_space-1);
                this.d_basis_polar   = DiscreteLP(this.d_number_polar-1);
                % normal LP expansion; we'll look at DPn later.
                this.d_basis_azimuth = DiscreteLP(2*this.d_number_azimuth-1);
                
            else
                error('Only DIM = 1, 2  supported')
            end
            
            
        end
        
        % ======================================================================
        %> @brief Update the boundary flux.
        % ======================================================================
        function this = update(this)
            
            % Boundary fluxes are stored as follows
            %   this.d_boundary_flux{side}(:, angle, this.d_g) 
            % i.e. for a given side, the data is stored by space and 
            % discrete angle and group.  The angle is the *cardinal*
            % index.  The quadrature is arranged (pol1,azi1,pol1,azi2,...)

            % Compute the incident flux if needed.
            if this.d_set == 0 && group(this.d_boundary) == this.d_order_group
                this.d_set = 1;   
            
                Ps = this.d_basis_space(:, this.d_order_space+1);
                Pa = this.d_basis_azimuth(:, this.d_order_azimuth+1);
                Pp = this.d_basis_polar(:, this.d_order_polar+1);
                
                for o = 1:length(this.d_octants(:, 1))
                    
                    o_in  = this.d_octants(o, 1); % incident octant
                    
                    f = zeros(this.d_number_space, ...
                        this.d_number_azimuth*this.d_number_polar);

                    if o == 1
                        a1 = 1; 
                        a2 = this.d_number_azimuth;
                        a3 = 1;
                    else
                        a1 = this.d_number_azimuth;
                        a2 = 1;
                        a3 = -1;
                    end
                    for s = 1:this.d_number_space
                        angle = 1;
                        for a = a1:a3:a2
                            for p = 1:this.d_number_polar
                                if 1 == 1%s == 15 && a == 1
                                f(s, angle) = Ps(s)*...
                                    Pa(a+(o-1)*this.d_number_azimuth)*Pp(p);
                                end
                                angle = angle + 1;
                            end
                        end
                    end        
                    
                    if this.d_side == Mesh.LEFT || this.d_side == Mesh.RIGHT
                        set_psi_v_octant(this.d_boundary, o_in, f, Boundary.IN);
                    elseif this.d_side == Mesh.TOP || this.d_side == Mesh.BOTTOM
                        set_psi_h_octant(this.d_boundary, o_in, f, Boundary.IN);
                    end
                    
                end
            end
            
        end
        
    end
    
end