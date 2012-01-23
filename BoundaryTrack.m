%> @file  BoundaryTrack.m
%> @brief BoundaryTrack class definition.
% ==============================================================================
%> @brief Holds the incident and outgoing boundary fluxes for tracks.
%
%> More here...
% ==============================================================================
classdef BoundaryTrack < Boundary
    
    methods
        
        % ======================================================================
        %> @brief Class constructor
        %>
        %> All boundary flux arrays are sized, and the boundary conditions are
        %> constructed for each surface.
        %>
        %> @param input         User input.
        %> @param track         Spatial track.
        %> @param quadrature  	Angular track.
        %>
        %> @return Instance of the BoundaryTrack class.
        % ======================================================================
        function this = BoundaryTrack(input, track, quadrature)
        
            % Call the base class.
            this = this@Boundary(input, track, quadrature);
            

            this.d_boundary_flux = cell(2*track.DIM, 1);
            
            
            ng = get(input, 'number_groups');
            
            
            this.d_boundary_flux{track.LEFT} = ...
                cell(number_azimuth(quadrature), 1);
            this.d_boundary_flux{track.RIGHT} = ...
                cell(number_azimuth(quadrature), 1);
            this.d_boundary_flux{track.TOP} = ...
                cell(number_azimuth(quadrature), 1);
            this.d_boundary_flux{track.BOTTOM} = ...
                cell(number_azimuth(quadrature), 1);
            
            for o = 1:4
                for m = 1:number_azimuth_octant(quadrature)
                    % Cardinal index.
                    mm = m + (o-1)*number_azimuth_octant(quadrature);
                    % LEFT BC's (phi, theta, v-surface track, group)
                    this.d_boundary_flux{track.LEFT}{mm} = ...
                        zeros(number_polar(quadrature), ...
                        number_y(quadrature, m), ng);
                    % RIGHT BC's (phi, theta, v-surface track, group)
                    this.d_boundary_flux{track.RIGHT}{mm} = ...
                        zeros(number_polar(quadrature),  ...
                        number_y(quadrature, m), ng);
                    % BOTTOM BC's (phi, theta, h-surface track, group)
                    this.d_boundary_flux{track.BOTTOM}{mm} = ...
                        zeros(number_polar(quadrature), ...
                        number_x(quadrature, m), ng);
                    % TOP BC's (phi, theta, h-surface track, group)
                    this.d_boundary_flux{track.TOP}{mm} = ...
                        zeros(number_polar(quadrature), ...
                        number_x(quadrature, m), ng);
                end
            end
        end
        
        % ======================================================================
        %> @brief Get a face incident/exiting boundary flux.
        %
        %> @param o     Octant index.
        %> @param a     Angle index (within an octant).
        %> @param inout Incoming/outgoing switch 
        %> @return      Boundary flux array.
        % ======================================================================
        function f = get_psi(this, o, a, side)
            angle = index(this.d_quadrature, o, a);
            f = this.d_boundary_flux{side}{angle}(:, :, this.d_g);
        end
        
        % ======================================================================
        %> @brief Set a vertical face incident boundary flux.
        %
        %> @param o     Octant index.
        %> @param a     Angle index (within an octant).
        %> @param f     Flux array.
        %> @param inout Incoming/outgoing switch         
        % ======================================================================
        function this = set_psi(this, o, a, side, f)
            angle = index(this.d_quadrature, o, a);
            this.d_boundary_flux{side}{angle}(:, :, this.d_g) = f;
        end
        
        
    end
    
end