%> @file  Response.m
%> @brief Response boundary condition class definition.
% ==============================================================================
%> @brief Isotropic boundary condition.
%
%> 
% ==============================================================================
classdef IsotropicBC < BoundaryCondition
    
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
        %> @return Instance of the IsotropicBC class.
        % ======================================================================
        function this = IsotropicBC(boundary, input, mesh, quadrature, side)
            this = this@BoundaryCondition(boundary, input, mesh, quadrature, side);

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
                end
               
            else
                error('Only DIM = 1, 2  supported')
            end
            

        end

        % ======================================================================
        %> @brief Initialize the boundary condition.
        %>
        %> For fixed boundaries, this is the only relevant call.
        % ======================================================================
        function this = initialize(this)   
            ng = get(this.d_input, 'number_groups');
            % Get the spectrum
            if (length(get(this.d_input, 'bc_isotropic_spectrum')) == ng)
                spectrum = get(this.d_input, 'bc_isotropic_spectrum');
            else
                spectrum = zeros(get(this.d_input, 'number_groups'), 1);
            end
            for g = 1:ng
                
                set_group(this.d_boundary, g);
                
                for o = 1:length(this.d_octants(:, 1))
                    
                    o_in  = this.d_octants(o, 1); % Put fluxes IN this octant
                    %o_out = this.d_octants(o, 2); % Get fluxes OUT of this octant
                    
                    if this.d_side == Mesh.LEFT || this.d_side == Mesh.RIGHT
                        
                        f = spectrum(g)*ones(size(get_psi_v_octant(this.d_boundary, o_in, Boundary.OUT)));
                        set_psi_v_octant(this.d_boundary, o_in, f, Boundary.IN);
                        
                    elseif this.d_side == Mesh.TOP || this.d_side == Mesh.BOTTOM
                        
                        f = spectrum(g)*ones(size(get_psi_h_octant(this.d_boundary, o_in, Boundary.OUT)));
                        set_psi_h_octant(this.d_boundary, o_in, f, Boundary.IN);
                        
                    else
                        error(' 3-D not done yet')
                    end
                    
                end
            end
            
        end
        
        % ======================================================================
        %> @brief Update the boundary flux.
        % ======================================================================
        function this = update(this)
            % Nothing to do.
            initialize(this);
        end
        
    end
    
end