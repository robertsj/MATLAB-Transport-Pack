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
        %> Angular expansion order.  Maximum is number angles - 1.
        d_order_angle = 14
        %> Spatial expansion order.  Maximum is number of points - 1.
        d_order_space = 27
        %> Discrete Legendre polynomial spatial basis
        d_basis_space
        %> Discrete Legendre polynomial angular basis
        d_basis_angle
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
        function this = AppxReflectiveTrack(boundary, input, mesh, quadrature, side)
            this = this@BoundaryConditionMOC(boundary, input, mesh, quadrature, side); 
            
            %this.d_index = boundary.d_side_index{side};
            this.d_number_tracks = length(boundary.d_side_index{side}(:,1));
            build_index(this);
            build_basis(this);
            build_group(this);
            
            % Set the expansion orders as set by the user in the collocated
            % quadrature.
            this.d_order_angle = quadrature.d_order_angle;
            this.d_order_space = quadrature.d_order_space;
            
        end

        % ======================================================================
        %> @brief Set the boundary flux.
        %
        %> Note, nothing needs to be done for reflective, since the initial
        %> guess is zero.
        % ======================================================================
        function this = set(this)
            % do nothing
        end      
        
        % ======================================================================
        %> @brief Update the boundary flux.
        % ======================================================================
        function this = update(this)
            
            % Loop through all the tracks I own, grap the flux of each
            % track's feeder, and update my own track flux.
            psi_exact = zeros(this.d_number_tracks, 1);
            for i = 1:this.d_number_tracks
                
                for p = 1:number_polar(this.d_quadrature)
                    
                    % Get the feeder flux.
                    psi_exact(i) = ...
                        get_single_psi(this.d_boundary, this.d_index(i, 4), ...
                        this.d_index(i, 5), p, this.d_index(i, 6), 2);
                    
                end
                
                
                
            end

            % Approximate the flux.
            psi_appx = expand(this, psi_exact);
            
            % Debug plots.
            %             if length(psi_exact)==14
            %             a1 = 2; a2 = 1; a3 = 4;
            %             x1 = 1; x2 = 2; x3 = 3;
            %             elseif length(psi_exact)==62
            %             a1 = 2; a2 = 5; a3 = 8;
            %             x1 = 1; x2 = 5; x3 = 8;
            %             else
            %             a1 = 3; a2 = 7;  a3= 8; a4 = 12;
            %             x1 = 3; x2 = 14; x3 = 25;
            %             end
            %
            %             for i = 1:27
            %                pp       = psi_exact(this.d_group_x{i});
            %                phi(i)   = mean(pp);
            %             end
            %
            %
            %
            %             % Angles at POINTS
            %
            %             pe_a1 = psi_exact(this.d_group_x{x1});
            %             pa_a1 = psi_appx(this.d_group_x{x1});
            %             pe_a2 = psi_exact(this.d_group_x{x2});
            %             pa_a2 = psi_appx(this.d_group_x{x2});
            %             pe_a3 = psi_exact(this.d_group_x{x3});
            %             pa_a3 = psi_appx(this.d_group_x{x3});
            %
            %             % Points of an ANGLE
            %
            %             pe_x1 = psi_exact(this.d_group_a{a1});
            %             pa_x1 = psi_appx(this.d_group_a{a1});
            %             pe_x2 = psi_exact(this.d_group_a{a2});
            %             pa_x2 = psi_appx(this.d_group_a{a2});
            %             pe_x3 = psi_exact(this.d_group_a{a3});
            %             pa_x3 = psi_appx(this.d_group_a{a3});
            %             pe_x4 = psi_exact(this.d_group_a{a4});
            %             pa_x4 = psi_appx(this.d_group_a{a4});
            %
            % %             disp('side='),this.d_side
            %             figure(1+(this.d_side-1)*2)
            % %             plot(psi_exact)
            % %             lala = find(psi_exact);
            % %             if length(lala)
            % %             disp([' side ', num2str(this.d_side),' has flux at track ', num2str(lala)])
            % %             end
            %             subplot(4,1,1)
            %             plot(1:length(pe_a1),pe_a1,'k',1:length(pe_a1),pa_a1,'g--o','LineWidth', 2), grid on
            %             title('angle spectrum at location 1')
            %             subplot(4,1,2)
            %             plot(1:length(pe_a2),pe_a2,'k',1:length(pe_a2),pa_a2,'g--o','LineWidth', 2), grid on
            %             title('angle spectrum at location 2')
            %             subplot(4,1,3)
            %             plot(1:length(pe_a3),pe_a3,'k',1:length(pe_a3),pa_a3,'g--o','LineWidth', 2), grid on
            %             title('angle spectrum at location 3')
            %             subplot(4,1,4)
            %             plot(1:27,phi,'k','LineWidth', 2), grid on
            %             title('boundary average scalar flux')
            %
            %             figure(2+(this.d_side-1)*2)
            %             subplot(4,1,1)
            %             plot(1:length(pe_x1),pe_x1,'k',1:length(pe_x1),pa_x1,'g--o','LineWidth', 2), grid on
            %             title('angle 1 as a function of x')
            %             subplot(4,1,2)
            %             plot(1:length(pe_x2),pe_x2,'k',1:length(pe_x2),pa_x2,'g--o','LineWidth', 2), grid on
            %             title('angle 2 as a function of x')
            %             subplot(4,1,3)
            %             plot(1:length(pe_x3),pe_x3,'k',1:length(pe_x3),pa_x3,'g--o','LineWidth', 2), grid on
            %             title('angle 3 as a function of x')
            %             subplot(4,1,4)
            %             plot(1:length(pe_x4),pe_x4,'k',1:length(pe_x4),pa_x4,'g--o','LineWidth', 2), grid on
            %             title('angle 4 as a function of x')

            % Set the flux.
            for i = 1:this.d_number_tracks
                for p = 1:number_polar(this.d_quadrature)
                    
                    set_single_psi(this.d_boundary, this.d_index(i, 1), ...
                        this.d_index(i, 2), p, this.d_index(i, 3), ...
                        psi_appx(i), 1);
                    
                end
            end
            
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
            
            % Define how many tracks at each spacial location or 
            % in a given angle we can have.
            if this.d_quadrature.d_number_azimuth == 3
                this.d_n_space = [1 3];
                this.d_n_angle = [4 6];
            elseif this.d_quadrature.d_number_azimuth == 5
                this.d_n_space = [1 3  9];
                this.d_n_angle = [6 8 10];
                
            else % 7
                this.d_n_space = [1  3  9 27];
                this.d_n_angle = [8 10 12 14];
            end
            
            % Build the normalized DLP's.   
            this.d_basis_space = cell(length(this.d_n_space), 1);
            for i = 1:length(this.d_n_space)
                this.d_basis_space{i} = DiscreteLP(this.d_n_space(i)-1);
            end
            
            this.d_basis_angle = cell(length(this.d_n_angle), 1);
            for i = 1:length(this.d_n_angle);
                this.d_basis_angle{i} = DiscreteLP(this.d_n_angle(i)-1);
            end
            
        end
        
        % ======================================================================
        %> @brief Build space and angle index groups.
        %
        %> This builds two sets of indices.  The first contains a vector of
        %> all tracks at a particular spatial point, ordered by increasing
        %> angle.  Hence, an angular basis can be used directly on the 
        %> values at that point.
        %> 
        %> The second is a list of all tracks for a given angle, ordered by
        %> increasing (local) spatial location.  Hence, the spatial basis
        %> can be used directly for this set of tracks.
        %> 
        % ======================================================================
        function this = build_group(this)
           
            if this.d_quadrature.d_number_azimuth == 3
                
                % Tracks at a given origin.
                this.d_group_x = cell(3, 1);
                this.d_group_x{1} = [  2 5  8 11   ];
                this.d_group_x{2} = [1 3 6  9 12 14];
                this.d_group_x{3} = [  4 7 10 13   ];
                
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
                base = [5 14 23 32 41 50];
                this.d_group_x{1} = [      base           ];
                this.d_group_x{2} = [   2  base+1   59    ];
                this.d_group_x{3} = [      base+2         ];
                this.d_group_x{4} = [      base+3         ];
                this.d_group_x{5} = [ 1 3  base+4   60 62 ];
                this.d_group_x{6} = [      base+5         ];
                this.d_group_x{7} = [      base+6         ];
                this.d_group_x{8} = [   4  base+7   61    ];
                this.d_group_x{9} = [      base+8         ];
                
                % Tracks in a given direction.
                this.d_group_a    = cell(10, 1);
                this.d_group_a{1} = [1];
                this.d_group_a{2} = [2 3 4];
                this.d_group_a{3} = [ 4 + [1:9]];
                this.d_group_a{4} = [13 + [1:9]];
                this.d_group_a{5} = [22 + [1:9]];
                this.d_group_a{6} = [31 + [1:9]];
                this.d_group_a{7} = [40 + [1:9]];
                this.d_group_a{8} = [49 + [1:9]];
                this.d_group_a{9} = [59 60 61];
                this.d_group_a{10}= [62];
                
            elseif  this.d_quadrature.d_number_azimuth == 7
                
                % Tracks at a given origin.
                this.d_group_x     = cell(27, 1);
                base = 14 + 27.*[0 1 2 3 4 5 6 7];
                
                this.d_group_x{ 1} = base;
                this.d_group_x{ 2} = [     5  base+1   230];        	%*
                this.d_group_x{ 3} = [base+2]; 
                this.d_group_x{ 4} = [base+3];
                this.d_group_x{ 5} = [   2 6  base+4   231 239 ];       %**
                this.d_group_x{ 6} = [base+5];
                this.d_group_x{ 7} = [base+6];
                this.d_group_x{ 8} = [     7  base+7   232];            %*
                this.d_group_x{ 9} = [base+8];
                this.d_group_x{10} = [base+9];
                this.d_group_x{11} = [     8  base+10  233 ];           %*
                this.d_group_x{12} = [base+11];
                this.d_group_x{13} = [base+12];
                this.d_group_x{14} = [ 1 3 9  base+13  234 240 242 ];	%***        
                this.d_group_x{15} = [base+14];
                this.d_group_x{16} = [base+15];
                this.d_group_x{17} = [     10 base+16  235];           	%*
                this.d_group_x{18} = [base+17];
                this.d_group_x{19} = [base+18];
                this.d_group_x{20} = [     11 base+19  236];        	%*
                this.d_group_x{21} = [base+20];
                this.d_group_x{22} = [base+21];
                this.d_group_x{23} = [   4 12 base+22  237 241];     	%**
                this.d_group_x{24} = [base+23];
                this.d_group_x{25} = [base+24];
                this.d_group_x{26} = [     13 base+25  238];           	%*
                this.d_group_x{27} = [base+26];
              
                % Tracks in a given direction.
                this.d_group_a    = cell(14, 1);
                i = 1;
                this.d_group_a{1}   = [i:i+0 ]; i = i+1;
                this.d_group_a{2}   = [i:i+2 ]; i = i+3;
                this.d_group_a{3}   = [i:i+8 ]; i = i+9;
                this.d_group_a{4}   = [i:i+26]; i = i+27;
                this.d_group_a{5}   = [i:i+26]; i = i+27;
                this.d_group_a{6}   = [i:i+26]; i = i+27;
                this.d_group_a{7}   = [i:i+26]; i = i+27;
                this.d_group_a{8}   = [i:i+26]; i = i+27;
                this.d_group_a{9}   = [i:i+26]; i = i+27;
                this.d_group_a{10}  = [i:i+26]; i = i+27;
                this.d_group_a{11}  = [i:i+26]; i = i+27;
                this.d_group_a{12}  = [i:i+8 ]; i = i+9;         
                this.d_group_a{13}  = [i:i+2 ]; i = i+3;
                this.d_group_a{14}  = [i:i+0 ]; 
                assert(i==242);
                
            end
            
            
        end
        
        % ======================================================================
        %> @brief Expand the reflected flux and approximate.
        %
        %> This takes the exact reflected flux and replaces it with an
        %> (n,m)th order approximation in space and angle.
        %> 
        %> @todo Study the effect of approximating space or angle first.
        %>       It doesn't matter in a full expansion, but it does for
        %>       low orders as they are not independent variables.
        %> 
        % ======================================================================
        function psi_appx = expand(this, psi_exact)
           
            psi_appx = 0*psi_exact;
            
            % Expand in ANGLE at each spatial point.
            for i = 1:length(this.d_group_x)
                % Psi indices
                idx = this.d_group_x{i};
                % Number of values
                num   = length(idx);
                % Get the basis.
                P = this.d_basis_angle{this.d_n_angle==num};
                f  = psi_exact(idx);
                % Expansion coefficients.
                fc = f' * P;
                f  = 0*f;
                % M-th order angular approximation.  If the requested order
                % is higher than the number of points, the maximum possible
                % is used instead. (We might have 14 possible angles on a
                % surface, but the most accute angle might not be at this
                % point)
                for o = 1:min(num, this.d_order_angle+1) % isotropic
                    f = f + fc(o) * P(:, o);
                end
                psi_appx(idx) = f;
            end
            
            % Expand in SPACE for each angle.
            for i = 1:length(this.d_group_a)
                % Psi indices
                idx = this.d_group_a{i};
                % Number of values
                num   = length(idx);
                % Get the basis.
                P = this.d_basis_space{this.d_n_space==num};
                f  = psi_appx(idx);
                % Expansion coefficients.
                fc = f' * P;
                f = 0*f;
                % N-th order spatial approximation.
                for o = 1:min(num, this.d_order_space+1) % isotropic
                    f = f + fc(o) * P(:, o);
                end
                psi_appx(idx) = f;
            end

        end

    end
 
end

