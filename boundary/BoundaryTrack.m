%> @file  BoundaryTrack.m
%> @brief BoundaryTrack class definition.
% ==============================================================================
%> @brief Holds the incident and outgoing boundary fluxes for tracks.
%
%> More here...
% ==============================================================================
classdef BoundaryTrack < handle
    
    properties (Constant)
        INCIDENT = 1;
        EXITING  = 2;
    end
    properties (Access = public)
        %> Spatial track.
        d_track
        %> Angular track.
        d_quadrature
        %> Boundary flux cell array
        d_boundary_flux
        %> Cell array of \ref BoundaryCondition thisects for each surface
        d_bc
        %> Current group.
        d_g = 0;
    end
    properties
        d_octants
        d_index
        d_bindex
        d_side_index
    end
        
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
            this.d_track         = track;
            this.d_quadrature    = quadrature;

            % Build useful indices.
            build_reflect_index(this);
            build_side_index(this);
            
            % Boundary conditions
            this.d_bc = cell(2*track.DIM);
            
            if strcmp(get(input, 'bc_left'), 'vacuum')
                this.d_bc{track.LEFT}   = ...
                    VacuumTrack(this, track, quadrature, track.LEFT);
            elseif strcmp(get(input, 'bc_left'), 'reflect')
                this.d_bc{track.LEFT}   = ...
                    ReflectiveTrack(this, track, quadrature, track.LEFT);
         	elseif strcmp(get(input, 'bc_left'), 'reflect_r')
                this.d_bc{track.LEFT}   = ...
                    AppxReflectiveTrack(this, track, quadrature, track.LEFT);
            end
            
            if strcmp(get(input, 'bc_right'), 'vacuum')
                this.d_bc{track.RIGHT}   = ...
                    VacuumTrack(this, track, quadrature, track.RIGHT);
            elseif strcmp(get(input, 'bc_right'), 'reflect')
                this.d_bc{track.RIGHT}   = ...
                    ReflectiveTrack(this, track, quadrature, track.RIGHT);
            elseif strcmp(get(input, 'bc_right'), 'reflect_r')
                this.d_bc{track.RIGHT}   = ...
                    AppxReflectiveTrack(this, track, quadrature, track.RIGHT);
            end
            
            if strcmp(get(input, 'bc_bottom'), 'vacuum')
                this.d_bc{track.BOTTOM}   = ...
                    VacuumTrack(this, track, quadrature, track.BOTTOM);
            elseif strcmp(get(input, 'bc_bottom'), 'reflect')
                this.d_bc{track.BOTTOM}   = ...
                    ReflectiveTrack(this, track, quadrature, track.BOTTOM);
            elseif strcmp(get(input, 'bc_bottom'), 'reflect_r')
                this.d_bc{track.BOTTOM}   = ...
                    AppxReflectiveTrack(this, track, quadrature, track.BOTTOM);
            end
            
            if strcmp(get(input, 'bc_top'), 'vacuum')
                this.d_bc{track.TOP}   = ...
                    VacuumTrack(this, track, quadrature, track.TOP);
            elseif strcmp(get(input, 'bc_top'), 'reflect')
                this.d_bc{track.TOP}   = ...
                    ReflectiveTrack(this, track, quadrature, track.TOP);
            elseif strcmp(get(input, 'bc_top'), 'reflect_r')
                this.d_bc{track.TOP}   = ...
                    AppxReflectiveTrack(this, track, quadrature, track.TOP);
            end
            
            % Number of groups
            ng = get(input, 'number_groups');

            % Inititalize boundary fluxes.
            this.d_boundary_flux = cell(2, 1);
            this.d_boundary_flux{1} = cell(number_angles(quadrature), 1);
            this.d_boundary_flux{2} = cell(number_angles(quadrature), 1);
            
            for o = 1:4
                for m = 1:number_angles_octant(quadrature)
                    
                    % Cardinal index.
                    mm = m + (o-1)*number_angles_octant(quadrature);
                    
                    % INCIDENT FLUX (phi, theta, v-surface track, group)
                    this.d_boundary_flux{1}{mm} = ...
                        zeros(number_polar(quadrature), ...
                              number_tracks(quadrature, m), ...
                              ng);
                    % EXITING FLUX (phi, theta, v-surface track, group)
                    this.d_boundary_flux{2}{mm} = ...
                        zeros(number_polar(quadrature), ...
                              number_tracks(quadrature, m), ...
                              ng);                          
                end
            end
            IN  = 1;
            OUT = 2;
            
            this.d_octants(track.LEFT, 1, 1) =  1; % I go into the right, and
            this.d_octants(track.LEFT, 1, 2) =  2; % grab what goes left.
            this.d_octants(track.LEFT, 2, 1) =  4; % That is, I go into (I, IV) and
            this.d_octants(track.LEFT, 2, 2) =  3; % get what goes out of (II, III)
            
            this.d_octants(track.RIGHT, 1, 1) =  2; % I go into the left, and
            this.d_octants(track.RIGHT, 1, 2) =  1; % grab what goes right.
            this.d_octants(track.RIGHT, 2, 1) =  3; % That is, I go into (II, III) and
            this.d_octants(track.RIGHT, 2, 2) =  4; % get what goes out of (I, IV)
            
            this.d_octants(track.BOTTOM, 1, 1) =  1; % I go to the top, which is
            this.d_octants(track.BOTTOM, 2, 1) =  2; % octants 1 and 2
            this.d_octants(track.BOTTOM, 1, 2) =  4; % and get reflection
            this.d_octants(track.BOTTOM, 2, 2) =  3; % from 4 and 3, respectively
            
            this.d_octants(track.TOP, 1, 1) =  3; % I go to the bottom, which is
            this.d_octants(track.TOP, 2, 1) =  4; % octants 3 and 4
            this.d_octants(track.TOP, 1, 2) =  2; % and get reflection
            this.d_octants(track.TOP, 2, 2) =  1; % from 2 and 1, respectively


            
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
        function this = initialize(this, g)
            this.d_g = g;
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
        function this = set(this)
            for side = 1:2*this.d_track.DIM;
                update(this.d_bc{side});
            end
        end   
        
        % ======================================================================
        %> @brief Get incident or exiting flux for an angle.
        %
        %> @param o     Octant index.
        %> @param a     Angle index (within an octant).
        %> @param inout Incoming/outgoing switch 
        %> @return      Boundary flux array (polar, track)
        % ======================================================================
        function f = get_psi(this, o, a, inout)
            angle = index(this.d_quadrature, o, a);
            f(:, :) = this.d_boundary_flux{inout}{angle}(:, :, this.d_g);
        end
        
        % ======================================================================
        %> @brief Set incident or exiting flux for an angle.
        %
        %> @param o     Octant index.
        %> @param a     Angle index (within an octant).
        %> @param f     Flux array (polar, track)
        %> @param inout Incoming/outgoing switch         
        % ======================================================================
        function this = set_psi(this, o, a, f, inout)
            angle = index(this.d_quadrature, o, a);
            this.d_boundary_flux{inout}{angle}(:, :, this.d_g) = f(:, :);
        end
        
        % ======================================================================
        %> @brief Get incident or exiting flux for an angle and track
        %
        %> @param o     Octant index.
        %> @param a     Angle index (within an octant).
        %> @param inout Incoming/outgoing switch 
        %> @return      Boundary flux array (polar, track)
        % ======================================================================
        function f = get_single_psi(this, o, a, p, t, inout)
            angle = index(this.d_quadrature, o, a);
            f = this.d_boundary_flux{inout}{angle}(p, t, this.d_g);
        end
        
        % ======================================================================
        %> @brief Set incident or exiting flux for an angle and track
        %
        %> @param o     Octant index.
        %> @param a     Angle index (within an octant).
        %> @param f     Flux array (polar, track)
        %> @param inout Incoming/outgoing switch         
        % ======================================================================
        function this = set_single_psi(this, o, a, t, p, f, inout)
            angle = index(this.d_quadrature, o, a);
            this.d_boundary_flux{inout}{angle}(t, p, this.d_g) = f;
        end
        
        % ======================================================================
        %> @brief Get incident or exiting flux for all angles
        %
        %> @param inout Incoming/outgoing switch 
        %> @return      Boundary flux array {azimuth}(polar, track)
        % ======================================================================
        function f = get_boundary_flux(this, inout)
            f = this.d_boundary_flux{inout};
        end
        
        % ======================================================================
        %> @brief Get incident or exiting flux for all angles for a side.
        %
        %> @param inout Incoming/outgoing switch 
        %> @return      Boundary flux array {azimuth}(polar, track)
        % ======================================================================
        function f = get_side_flux(this, side, polar, inout) 
            for i = 1:2
                o = this.d_octants(side, i, inout);
                for a = 1:number_angles_octant(quadrature)
                    angle = index(this.d_quadrature, o, a);
                    a2 = a + o*number_angles_octant(quadrature);
                    if ((side < 3 && inout == 1) || (side > 2 && inout == 2))
                        for t = 1:number_y(quadrature, a);
                            f{a2}(t) = this.d_boundary_flux{inout}{angle};
                        end
                    end
                end
            end
        end
        
        % ======================================================================
        %> @brief Set incident or exiting flux for all angles for a side.
        %
        %> @param inout Incoming/outgoing switch 
        %> @return      Boundary flux array {azimuth}(polar, track)
        % ======================================================================
        function this = set_side_flux(this, side, f, inout)     
           this.d_boundary_flux{inout} = f;
        end 
        
        function bc = get_bc(this, side)
            bc = this.d_bc{side};
        end
        
        % ======================================================================
        %> @brief Build index of (o, a) -> (o', a') pairs and its inverse.
        % ======================================================================
        function this = build_reflect_index(this)
            index = cell(4, number_angles_octant(this.d_quadrature));
            % mirrored octants
            oct = [4 2; 3 1; 2 4; 1 3];
            for o = 1:4
                for a = 1:number_angles_octant(this.d_quadrature)
                    nx = number_x(this.d_quadrature, a);
                    ny = number_y(this.d_quadrature, a);
                    n  = nx + ny;
                    index{o, a} = zeros(n, 3);
                    for t = 1:nx
                        index{o, a}(t, 1) = oct(o, 1);
                        index{o, a}(t, 2) = a;
                        index{o, a}(t, 3) = ny + t;
                    end
                    for t = 1:ny
                        index{o, a}(t+nx, 1) = oct(o, 2);
                        index{o, a}(t+nx, 2) = a;
                        index{o, a}(t+nx, 3) = t; 
                    end
                end
            end
            bindex = index;
            for o = 1:4
                for a = 1:number_angles_octant(this.d_quadrature)
                    nx = number_x(this.d_quadrature, a);
                    ny = number_y(this.d_quadrature, a);
                    n  = nx + ny;
                    for t = 1:n
                        oo = index{o, a}(t, 1);
                        aa = index{o, a}(t, 2);
                        tt = index{o, a}(t, 3);
                        bindex{oo, aa}(tt, 1) = o;
                        bindex{oo, aa}(tt, 2) = a;
                        bindex{oo, aa}(tt, 3) = t;
                    end
                end
            end
            this.d_index = index;
            this.d_bindex = bindex;
        end
        
        % ======================================================================
        %> @brief Build index of (o, a, t) for each side.
        %
        %> For a given side, have an array of three columns listing o, a,
        %> and t that it has as input.  Note, for each side, the tracks are
        %> listed by angle-space, from (pi, 0) and (i=1, i=# points).  This
        %> facilitates defining functions on a boundary as a function of
        %> angle or space.
        % ======================================================================
        function this = build_side_index(this)
            % incident octants
            oct = [1 4; 3 2; 2 1; 4 3];  % left, right from side perspective
            
            side_index = cell(4, 1);
            for s = 1:2
                side_index{s} = zeros(1, 3);
                count = 1;
                for o = 1:2
                    oi = oct(s, o);
                    a1 = 1; a2 = number_angles_octant(this.d_quadrature); a3=1;
                    if o == 1
                        a1 = a2; a2 = 1; a3 = -1; 
                    end
                    for a = a1:a3:a2
                        ny = number_y(this.d_quadrature, a);
                        b1 = 1; b2 = ny; b3 = 1;
                        if o == 2
                            b1 = b2; b2 = 1; b3 = -1;
                        end
                        for t = b1:b3:b2
                            side_index{s}(count, 1) = oi;
                            side_index{s}(count, 2) = a;
                            side_index{s}(count, 3) = t;
                            count = count + 1;
                        end
                    end
                end
            end
            for s = 3:4
                side_index{s} = zeros(1, 3);
                count = 1;
                for o = 1:2
                    oi = oct(s, o);
                    a1 = 1; a2 = number_angles_octant(this.d_quadrature); a3=1;
                    if o == 2
                        a1 = a2; a2 = 1; a3 = -1; 
                    end
                    for a = a1:a3:a2
                        nx = number_x(this.d_quadrature, a);
                        ny = number_y(this.d_quadrature, a);
                        n  = nx + ny;
                        b1 = ny+1; b2 = n; b3 = 1;
                        if o == 1
                            b2 = b1; b1 = n; b3 = -1;
                        end
                        for t = b1:b3:b2
                            side_index{s}(count, 1) = oi;
                            side_index{s}(count, 2) = a;
                            side_index{s}(count, 3) = t;
                            count = count + 1;
                        end
                    end
                end
            end
            this.d_side_index = side_index;
        end
    end
    
    
end