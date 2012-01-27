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
        %> Index of all tracks as a function of local angle index (increasing x)
        d_group_a
        %> Index of all tracks as a function of spatial point (increasing phi)
        d_group_x
        d_n_space
        d_n_angle
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
            ao = 1;  a = this.d_index(1, 5); tt = 0; aa=1;
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
        
        function this = build_group(this)
           
            if this.d_quadrature.d_number_azimuth == 3
                
                % Tracks at a given origin.
                this.d_group_x = cell(3, 1);
                this.d_group_x{1} = [2 5 8 11];
                this.d_group_x{2} = [1 3 6 9 12 14];
                this.d_group_x{3} = [4 7 10 13];
                
                % Tracks in a given direction.
                this.d_group_a = cell(6, 1);
                this.d_group_a{1} = [1];
                this.d_group_a{2} = [2 3 4];
                this.d_group_a{3} = [5 6 7];
                this.d_group_a{4} = [8 9 10];
                this.d_group_a{5} = [11 12 13];
                this.d_group_a{6} = [14];
                
            elseif this.d_quadrature.d_number_azimuth == 5
                
                % Tracks at a given origin.
                this.d_group_x    = cell(9, 1);
                this.d_group_x{1} = [5 14 23 32 41 50];
                this.d_group_x{2} = [2 this.d_group_x{1}+1 59];
                this.d_group_x{3} = [this.d_group_x{1}+2];
                this.d_group_x{4} = [this.d_group_x{1}+3];
                this.d_group_x{5} = [1 3 this.d_group_x{1}+4 60 62 ];
                this.d_group_x{6} = [this.d_group_x{1}+5];
                this.d_group_x{7} = [this.d_group_x{1}+6];
                this.d_group_x{8} = [4 this.d_group_x{1}+7 61];
                this.d_group_x{9} = [this.d_group_x{1}+8];
                
                % Tracks in a given direction.
                this.d_group_a    = cell(10, 1);
                this.d_group_a{1} = [1];
                this.d_group_a{2} = [2 3 4];
                this.d_group_a{3} = [ 4 + 1:9];
                this.d_group_a{4} = [13 + 1:9];
                this.d_group_a{5} = [22 + 1:9];
                this.d_group_a{6} = [31 + 1:9];
                this.d_group_a{7} = [40 + 1:9];
                this.d_group_a{8} = [49 + 1:9];
                this.d_group_a{9} = [59 60 61];
                this.d_group_a{10}= [62];
                
            elseif  this.d_quadrature.d_number_azimuth == 7
                
                % Tracks at a given origin.
                this.d_group_x     = cell(27, 1);
                this.d_group_x{ 1} = [5 14 23 32 41 50];
                this.d_group_x{ 2} = [2 this.d_group_x{1}+1 59];
                this.d_group_x{ 3} = [this.d_group_x{1}+2];
                this.d_group_x{ 4} = [this.d_group_x{1}+3];
                this.d_group_x{ 5} = [1 3 this.d_group_x{1}+4 60 62 ];
                this.d_group_x{ 6} = [this.d_group_x{1}+5];
                this.d_group_x{ 7} = [this.d_group_x{1}+6];
                this.d_group_x{ 8} = [4 this.d_group_x{1}+7 61];
                this.d_group_x{ 9} = [this.d_group_x{1}+8];
                this.d_group_x{10} = [5 14 23 32 41 50];
                this.d_group_x{11} = [2 this.d_group_x{1}+1 59];
                this.d_group_x{12} = [this.d_group_x{1}+2];
                this.d_group_x{13} = [this.d_group_x{1}+3];
                this.d_group_x{14} = [this.d_group_x{1}+8];                
                this.d_group_x{15} = [1 3 this.d_group_x{1}+4 60 62 ];
                this.d_group_x{16} = [this.d_group_x{1}+5];
                this.d_group_x{17} = [this.d_group_x{1}+6];
                this.d_group_x{18} = [4 this.d_group_x{1}+7 61];
                this.d_group_x{19} = [this.d_group_x{1}+8];
                this.d_group_x{20} = [5 14 23 32 41 50];
                this.d_group_x{21} = [this.d_group_x{1}+8];
                this.d_group_x{22} = [2 this.d_group_x{1}+1 59];
                this.d_group_x{23} = [this.d_group_x{1}+2];
                this.d_group_x{24} = [this.d_group_x{1}+3];
                this.d_group_x{25} = [1 3 this.d_group_x{1}+4 60 62 ];
                this.d_group_x{26} = [this.d_group_x{1}+5];
                this.d_group_x{27} = [this.d_group_x{1}+6];
              
                % Tracks in a given direction.
                this.d_group_a    = cell(14, 1);
                this.d_group_a{1} = [1];
                this.d_group_a{2} = [2 3 4];
                this.d_group_a{3} = [ 4 + 1:9];
                this.d_group_a{4} = [13 + 1:9];
                this.d_group_a{5} = [22 + 1:9];
                this.d_group_a{6} = [31 + 1:9];
                this.d_group_a{7} = [40 + 1:9];
                this.d_group_a{8} = [49 + 1:9];
                this.d_group_a{9} = [40 + 1:9];
                this.d_group_a{10} = [49 + 1:9];
                this.d_group_a{11} = [40 + 1:9];
                this.d_group_a{12} = [49 + 1:9];                
                this.d_group_a{13} = [59 60 61];
                this.d_group_a{14}= [62];
                
            end
            
            
        end
        
        
        function psi_appx = expand(this, psi_exact)
           
            
        end
        
        
    end
 
end

