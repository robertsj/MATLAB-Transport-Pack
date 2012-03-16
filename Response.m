%> @file  Response.m
%> @brief Response boundary condition class definition.
% ==============================================================================
%> @brief Response boundary condition.
%
%> Currently, this is limited to 2-d problems.
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
                this.d_number_polar   = number_polar(this.d_quadrature);
                this.d_number_azimuth = number_azimuth(this.d_quadrature);                
                
                this.d_basis_space   = DiscreteLP(this.d_number_space-1);
                this.d_basis_polar   = DiscreteLP(this.d_number_polar-1);
                % normal LP expansion; we'll look at DPn later.
                this.d_basis_azimuth = DiscreteLP(2*this.d_number_azimuth-1);
                
                % This azimuth basis is for left-to-right.  The angles as stored
                % are in mu_1 < mu_2 < ... < mu_N (or eta_1 > ...)
                % Hence, they are in order for vertical surfaces.
                if side > 2
                    b  = this.d_basis_azimuth;
                    lb = length(b);
                    this.d_basis_azimuth(1:lb/2, :) = b(lb/2:-1:1, :);
                    this.d_basis_azimuth(lb/2+1:end, :) = b(end:-1:lb/2+1, :);
                end
            else
                error('Only DIM = 1, 2  supported')
            end
            
            
        end
        
        % ======================================================================
        %> @brief Set the boundary flux.
        %
        %> Since the response condition is fixed, we set it once.
        % ======================================================================
        function this = set(this)

            % Boundary fluxes are stored as follows
            %   this.d_boundary_flux{side}(:, angle, this.d_g) 
            % i.e. for a given side, the data is stored by space and 
            % discrete angle and group.  The angle is the *cardinal*
            % index.  The quadrature is arranged (pol1,azi1,pol1,azi2,...)

            % Set the group of the incident response.
            set_group(this.d_boundary, this.d_order_group);
          
            % Set incident bases.
            Ps = this.d_basis_space(:, this.d_order_space+1);
            Pa = this.d_basis_azimuth(:, this.d_order_azimuth+1);
            Pp = this.d_basis_polar(:, this.d_order_polar+1);
                
            for o = 1:length(this.d_octants(:, 1))
                
                o_in  = this.d_octants(o, 1); % incident octant
                
                f = zeros(this.d_number_space, ...
                    this.d_number_azimuth*this.d_number_polar);
                    
                % If oo = 1, it's a vertical side, else a horizontal side.
                % The indices are selected so that psi is defined as a
                % function of monotonically increasing angle with respect to
                % its incident surface.  In other words, all incident
                % angular fluxes are defined left to right along an edge
                % both in space and angle.
                
                if (o == 1) % this makes angles symmetric
                    a1 = 1;
                    a2 = this.d_number_azimuth;
                    a3 = 1;
                else
                    a1 = this.d_number_azimuth;
                    a2 = 1;
                    a3 = -1;
                end

                
                for s = 1:this.d_number_space
                    ang = 1; % cardinal index
                    for a = a1:a3:a2
                        for p = 1:this.d_number_polar
                            f(s, (a-1)*this.d_number_polar+p) =  Ps(s)*...
                                Pa(ang+(o-1)*this.d_number_azimuth)*Pp(p);
                            
                        end
                        ang = ang + 1;
                    end
                end
                
                if this.d_side == Mesh.LEFT || this.d_side == Mesh.RIGHT
                    set_psi_v_octant(this.d_boundary, o_in, f, Boundary.IN);
                elseif this.d_side == Mesh.TOP || this.d_side == Mesh.BOTTOM
                    set_psi_h_octant(this.d_boundary, o_in, f, Boundary.IN);
                end
                
            end
            
        end
        
        % ======================================================================
        %> @brief Update the boundary flux.
        %
        %> Note, nothing needs to be done, since we've already fixed the
        %> incident response.
        % ======================================================================
        function this = update(this)
            % do nothing
        end 
        
        % ======================================================================
        %> @brief (Re)set the incident orders.
        % ======================================================================
        function this = set_orders(this, go, so, ao, po)
             this.d_order_group   = go;
             this.d_order_space   = so;
             this.d_order_azimuth = ao;
             this.d_order_polar   = po;
        end             
        
    end
    
end