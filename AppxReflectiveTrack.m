%> @file  ReflectiveTrack.m
%> @brief Boundary condition class definition.
% ==============================================================================
%> @brief Approximate reflective MOC boundary condition class.
%
%> The reflective condition is approximated via an (n,m)th order expansion
%> in space(n) and azimuthal angle(m).
%>
% ==============================================================================
classdef AppxReflectiveTrack < BoundaryConditionMOC
    
    properties 
        d_index
        d_number_tracks
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
        %> @return Instance of the AppxReflectiveTrack class.
        % ======================================================================
        function this = AppxReflectiveTrack(boundary, mesh, quadrature, side)
            this = this@BoundaryConditionMOC(boundary, mesh, quadrature, side); 
            
            %this.d_index = boundary.d_side_index{side};
            this.d_number_tracks = length(boundary.d_side_index{side}(:,1));
            build_index(this);
            build_basis(this);
        end

        % ======================================================================
        %> @brief Update the boundary flux.
        % ======================================================================
        function this = update(this)
            oct = [1 4; 3 2; 2 1; 4 3];
            % Loop through all the tracks I own, grap the flux of each
            % track's feeder, and update my own track flux.
            ao = 1;  a = this.d_index(i, 5); tt = 0; aa=1;
            for i = 1:this.d_number_tracks
                
                for p = 1:number_polar(this.d_quadrature)
                    
                    tt = tt + 1;

                    
                    % Get the feeder flux.
                    psi = get_single_psi(this.d_boundary, this.d_index(i, 1), ...
                        this.d_index(i, 2), p, this.d_index(i, 3), 2);
                    
                    ao = a;
                    
                    
                    % my octant, angle, and track
                    a = this.d_index(i, 5);
                    if ao ~= a
                        tt = 1; aa = aa+1;
                    end
                    
                    psi_exact{aa}(tt, 1) = psi;
                    psi_exact{aa}(tt, 2) = tt;
                    
                    % Set my flux.
%                     set_single_psi(this.d_boundary, this.d_index(i, 4), ...
%                         this.d_index(i, 5), p, this.d_index(i, 6), psi, 1);
                
                end
                
            end

            psi_appx = expand(this, psi_exact);
        end
        
        
        
    end
    
    methods (Access = protected)
       
        % ======================================================================
        %> @brief Build simple table of what I grab and what I set.
        % ======================================================================
        function this = build_index(this)
            
            % Presize my index.
            this.d_index = zeros(this.d_number_tracks, 6);
            
            % Go through all tracks that leave me.
            for i = 1:this.d_number_tracks
                
                % These specify my tracks.
                o = this.d_boundary.d_side_index{this.d_side}(i, 1);
                a = this.d_boundary.d_side_index{this.d_side}(i, 2);
                t = this.d_boundary.d_side_index{this.d_side}(i, 3);
                
                % Get who feeds me.
                o_f = this.d_boundary.d_bindex{o, a}(t, 1);
                a_f = this.d_boundary.d_bindex{o, a}(t, 2);
                t_f = this.d_boundary.d_bindex{o, a}(t, 3);
                
                % Add to the table.
                this.d_index(i, :) = [o a t o_f a_f t_f];
            
                
            end
            
            
        end
        
        % ======================================================================
        %> @brief Build my basis.
        %
        %> For m azimuthal angles per quadrant, we need 1+log_3 m spatial
        %> basis sets.  For n=9 and n=27, we need 3 and 4 sets.  The sets
        %> have 1, 3, 9, and 27 spatial points.  The lone point isn't
        %> really a basis, though...
        %>
        %> We need 3 and 4 angle sets, too.  For n=27, there are 7 angles 
        %> per quadrant, or 14 over the relevant half-space.  Then 10, 6,
        %> 4, and 2.  
        %> 
        % ======================================================================
        function this = build_basis(this)
            
            if this.d_quadrature.d_number_azimuth == 3
                this.d_n_space = [1 3];
                this.d_n_angle = [4 6];
                A_i     = [1 2 1];
            elseif this.d_quadrature.d_number_azimuth == 9
                this.d_n_space = [1 3 9];
                this.d_n_angle = [6 8 10];
                
            else
                this.d_n_space = [1 3 9 27];
                this.d_n_angle = [8 10 12 14];
            end
            
            % Build the normalized DLP's.
            
            this.d_basis_space = cell(this.d_n_space, 1);
            for i = 1:length(this.d_n_space)
                this.d_basis_space{i} = DiscreteLP(this.d_n_space(i));
            end
            
            this.d_basis_angle = cell(this.d_n_angle, 1);
            for i = 1:length(this.d_n_angle);
                this.d_basis_space{i} = DiscreteLP(this.d_n_angle(i));
            end
            
        end
        
        
        function psi_appx = expand(this, psi_exact)
           
            psi_exact_2 = cell(this.d_quadrature.d_number_space, 1);
            for a = 1:length(psi_exact)
                nx = length(psi_exact{a});
                if nx == 1
                    x_basis_index = length(this.d_basis_space);
                    psi_exact{a}(:, 2) = psi_exact{a}(:, 2)*7 - 1;
                elseif nx == 3
                    x_basis_index = length(this.d_basis_space) - 1;
                    psi_exact{a}(:, 2) = psi_exact{a}(:, 2)*5 - 1;
                elseif nx == 9
                    x_basis_index = length(this.d_basis_space) - 2;
                    psi_exact{a}(:, 2) = psi_exact{a}(:, 2)*3 - 1;
                else
                    x_basis_index = length(this.d_basis_space) - 3;
                end

                % the exact vector for this angle over all space
                f   = psi_exact{a}(1, :);
                Px  = this.d_basis_space{x_basis_index};
                % the legendre moments in SPACE
                fx  = f' * Px;
                % the APPROXIMATE vector in SPACE for this angle
                ff  = f*0;
                
                n = min(nx, 2);
                for i = 1:n
                    ff = ff + Px(:, i)*fx(i);
                end
                for i = 1:nx
                    l = length(psi_exact_2{psi_exact{a}(i, 2)});
                    psi_exact_2{psi_exact{a}(i, 2)}(l+1) = ff(i);
                end                
            end
             
            
            for t = this.d_quadrature.d_number_space
                na = length(psi_exact_2{t});
                a_basis_index = find( this.d_n_angle == na);
                
                
            end
                
                
            end
            
        end
        
        
    end
 
end

