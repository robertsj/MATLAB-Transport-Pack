%> @file  PinCell.m
%> @brief PinCell class definition.
% ==============================================================================
%> @brief  Represents reactor assembly pin in two dimensions.
%
%> This class represents a pin cell, including its geometry and associated
%> materials definition.  Currently, only a single square coolant region
%> is allowed with a possible central pin with one or more concentric regions.
%> 
%> The PinCell class can be used as a mesh for discrete ordinates calculations
%> via the \ref meshify routine, as PinCell inherits from Mesh2D.
%> Additionally, PinCell has a built in tracking routine for solution by
%> the method of characteristics.  The tracking is limited to four way
%> rotation symmetry.
%>
% ==============================================================================
classdef PinCell < Mesh2D
    
    properties (Access = public)
        
        %> @name Basic Data
        %> @{
        %
        %> Pin cell pitch (i.e. width)
        d_pitch
        %> Center (local)
        d_center
        %> Pin radii
        d_radii
        %> Number of radii
        d_number_radii
        %> Number of regions
        d_number_regions
        %> Region material map
        d_region_mat_map
        %> Region volumes
        d_region_volume
        %
        %> @}
        
        %> @name Mesh Data
        %> @{
        %
        %> Am I meshed?
        %d_meshed = 0;
        %
        %> @}
        
        %> @name Track Data
        %> @{
        %
        %> Am I tracked>
        d_tracked
        %> Quadrature
        d_quadrature
        %> Entrance points for each angle
        d_enter
        %> Exit points for each angle
        d_exit
        %
        d_space
        %> Track widths.
        d_track_width
        %> Segment lengths, in order of tracking (and sweep)
        d_segment_length
        %> Segment region index.
        d_segment_region
        %> Segment exponential coefficients. (Not used yet).
        d_segment_coef
        %> Number of segments by (angle, track).  For debugging.
        d_number_segments
        %
        %> @}
          
    end
    
    methods (Access = public)
        
        % ======================================================================
        %> @brief Class constructor
        %
        %> The user specifies the box size and the radii of pins.  An
        %> empty radius vector ("[]")
        %> indicates a homogeneous cell.  The radii 
        %> must be given from smallest to largest.
        %> A vector of
        %> material identifiers must also be given.  The size of this must
        %> be equal to the number of regions defined by the box and radii.
        %> The are arranged in-to-out. 
        %>
        %> Finally, the user can subdivide the domain.  This a level
        %> of 0 does nothing.  A level of 1 divides the domain along the
        %> the x and y axes (though the origin.  A level of 2 adds
        %> dividing planes at +/- 45 degrees through the origin.
        %>
        %> @param  pitch 	Pin cell width
        %> @param  radii 	Fuel pin radii, if any.
        %> @param  mat_map	Material identifiers for each region.
        %> @param  subdiv   Level of subcell division.
        %> @return Instance of the PinCell class.
        % ======================================================================
        function obj = PinCell(pitch, radii, mat_map, subdiv)
            % Precondition
            DBC.Require('length(radii) + 1 == length(mat_map)');
            
            obj.d_pitch          = pitch;
            obj.d_center         = [0.5*pitch 0.5*pitch];
            obj.d_radii          = radii;
            obj.d_number_radii   = length(radii);
            obj.d_number_regions = length(radii) + 1;
%            obj.d_sub_divide     = subdiv;
            
            % Don't allow radii >= to the pitch for simplicity.
            if max(radii) >= pitch/2
                error('Radii must be smaller than half pitch')
            end
            if sum(radii(1:end-1)>radii(2:end)) > 0
               error('non monotonic increasing radii') 
            end
            
            % Compute the volumes.
            obj.d_region_volume = zeros(obj.d_number_regions, 1);
            if obj.d_number_radii > 0
            	obj.d_region_volume(1) = pi*radii(1)^2;
            end
            for i = 2:obj.d_number_radii
                obj.d_region_volume(i) = ...
                    pi*radii(i)^2 - obj.d_region_volume(i-1);
            end
            if obj.d_number_radii > 0 
                obj.d_region_volume(end) = ...
                    pitch^2 - pi*radii(end)^2;
            else
                obj.d_region_volume(end) = pitch^2;
            end
            v1 = sum(obj.d_region_volume); 
            v2 = pitch^2;
            DBC.Ensure('v1 == v2');
            obj.d_region_mat_map = mat_map;  
            % I'm not meshed yet.
            obj.d_meshed = 0;
        end
        
        % ======================================================================
        %> @brief Construct a mesh in the pin cell.
        %
        %> This is a crude approach that given a number of meshes, a
        %> regular grid is defined.  Mesh materials are defined as the
        %> material in the mid point of a fine mesh. An alternative would
        %> be to define new materials via volume homogenization.
        %>
        %> @param  number_mesh  Number of uniformly spaces meshes per axis.
        % ======================================================================
        function obj = meshify(obj, number_meshes)
            
            width = obj.d_pitch / number_meshes;
            obj.d_dx  = width * ones(number_meshes, 1);
            obj.d_dy  = width * ones(number_meshes, 1);
            obj.d_number_cells_x = number_meshes;
            obj.d_number_cells_y = number_meshes;
            obj.d_number_cells   = number_meshes^2;
            
            tmp_material_map = zeros(number_meshes); 
            tmp_region_map   = zeros(number_meshes);
            
            volume = zeros(obj.d_number_regions, 1);
            for j = 1:number_meshes
                for i = 1:number_meshes
                    % Which region am I in?
                    r = find_region(obj, i, j);
                    % Which material does this region have?
                    m = obj.d_region_mat_map(r);
                    % Assign the values.
                    tmp_region_map(i, j) = r;
                    tmp_material_map(i, j)    = m;
                    % Add volume
                    volume(r) = volume(r) + width^2;
                end
            end
            
            % Add maps
            obj.d_mesh_map = containers.Map('MATERIAL', tmp_material_map);
            obj.d_mesh_map('REGION') = tmp_region_map;
            
            % Compute the relative error in volume.
            error_volume = (volume - obj.d_region_volume)./obj.d_region_volume;
            disp('Volume relative errors using mesh-center selection')
            for i = 1:length(volume)
                fprintf(' Reg: %2i, Vol. Rel. Err.: %12.8f\n', ...
                    i, error_volume(i));
            end
            obj.d_meshed = 1;
        end
        
        % ======================================================================
        %> @brief Perform tracking over the pin cell.
        %
        %> This is a crude approach that given a number of meshes, a
        %> regular grid is defined.  Mesh materials are defined as the
        %> material in the mid point of a fine mesh. An alternative would
        %> be to define new materials via volume homogenization.
        %>
        %> @param  q  MOC quadrature
        % ======================================================================
        function track(obj, q, mat)
            
            if obj.d_meshed
                error('Not doing meshing and tracking yet')
            else
                obj.d_tracked = 1;
            end
            obj.d_number_cells = obj.d_number_regions;
            pitch = obj.d_pitch;
            obj.d_quadrature = q;
            obj.d_enter = cell(number_azimuth(q), 1);
            obj.d_exit  = cell(number_azimuth(q), 1);
            
            for m = 1:number_angles_octant(q)  
                obj.d_enter{m} = q.d_enter{m}*pitch;
                obj.d_exit{m}  = q.d_exit{m}*pitch;
                obj.d_space(m) = q.d_space(m)*pitch;
            end
            
            % Tracking.  We make two passes.  The first pass counts the
            % the total number of tracks total.  With this, we size the
            % track length vector and region id vector.
            approximate_volume = zeros(obj.d_number_regions, 1);
            s = 1; % segment index
            
            
            for z = 1:number_angles_octant(q) 
                cos_phi = cos(phi(q, 1, z));
                tan_phi = tan(phi(q, 1, z));
                for t = 1:number_tracks(q, z)

                    % All points are translated by a half pitch to put
                    %  the cell center at the origin.
                    x_i = obj.d_enter{z}(t, 1)-pitch/2; % incident
                    y_i = obj.d_enter{z}(t, 2)-pitch/2; %   x and y
                    x_o = obj.d_exit{z}(t, 1)-pitch/2;  % outgoing
                    y_o = obj.d_exit{z}(t, 2)-pitch/2;  %   x and y
                    % total length
                    l   = sqrt( (x_o-x_i)^2 + (y_o-y_i)^2 ); 
                    % y intercept
                    b = y_i - (x_i * tan_phi);
                    % the closest approach to origin is
                    d_o = abs(cos_phi * b);
                    % how many cylinder regions to we intercept?
                    num_cyl = sum(obj.d_radii > d_o);
                    % number of segments 
                    num_seg = 2*num_cyl + 1;
                    
                    % Set segment region things.
                    obj.d_number_segments(z, t)=num_seg;
                    obj.d_segment_length{z}{t}=zeros(num_seg,1);
                    obj.d_segment_region{z}{t}=zeros(num_seg,1);
                    
                    ss = 1;
                    x = x_i;
                    y = y_i;
                    reg    = zeros(num_seg, 1);
                    reg(1) = obj.d_number_regions;
                    reg(end) = obj.d_number_regions;
                    if (num_seg > 1)
                        for i = 1:floor(num_seg/2)
                            idx = obj.d_number_radii - num_cyl + i;
                            reg(i+floor(num_seg/2))  = idx;
                            reg(ceil(num_seg/2)-i+1) = idx;
                        end
                        xx = zeros(2*num_cyl, 1); yy = xx;
                        for i = 1:num_cyl
                            % get outermost radius
                            r = obj.d_radii(end-i+1);
                            A = -tan_phi * b;
                            B = sqrt( tan_phi^2 * r^2 + r^2 - b^2);
                            C = tan_phi^2 + 1;
                            xx(i) = (A-B)/C;
                            yy(i) = tan_phi*xx(i) + b;
                            xx(2*num_cyl-i+1) = (A+B)/C;
                            yy(2*num_cyl-i+1) = tan_phi*xx(2*num_cyl-i+1)  + b;
                        end
                        x(ss+1:2*num_cyl+ss) = xx(:);
                        y(ss+1:2*num_cyl+ss) = yy(:);
                        
                    end
                    ss = ss+2*num_cyl+1;
                    x(ss) = x_o;
                    y(ss) = y_o;
                    
                    segl=zeros(num_seg, 1);
                    
                    for i = 1:num_seg
                        segl(i) = sqrt((x(i+1)-x(i))^2 + (y(i+1)-y(i))^2);
                        obj.d_segment_length{z}{t}(i) = segl(i);
                        obj.d_segment_region{z}{t}(i)  = reg(i);
                         
                        approximate_volume(reg(i)) = ...
                            approximate_volume(reg(i)) + ...
                            segl(i) * obj.d_space(z) * weight_phi(obj.d_quadrature, z) ;
%                        fprintf(' reg: %3i, segl: %12.8f, space: %12.8f , w: %12.8f, dvol  %12.8f \n', reg(i), segl(i) , obj.d_space(m) , weight_phi(obj.d_quadrature, m),  segl(i) * obj.d_space(m) * weight_phi(obj.d_quadrature, m));
%                         for j = 1:number_polar(q)
%                         	obj.d_segment_coef(s,j) = ...
%                                 exp(
%                         end
                        s = s + 1;
                    end
                    if ( abs(sum(segl)-l) > 1e-14 )
                        fprintf('%3i %3i %12.8f, %12.8f\n', z, t, sum(segl), l);
                        error(' error in track length ')
                    end
                end
                aaa=1;
            end
            approximate_volume = approximate_volume * 2/pi;
            % Fixup track lengths
            for z = 1:number_angles_octant(q) 
                for t = 1:number_tracks(q, z)
                    for i = 1:obj.d_number_segments(z, t)
                        region = obj.d_segment_region{z}{t}(i);
                        true_vol = obj.d_region_volume(region);
                        appx_vol = approximate_volume(region);
                        obj.d_segment_length{z}{t}(i) = ...
                            obj.d_segment_length{z}{t}(i) * true_vol/appx_vol;
                        
                    end
                end
            end
%             
%             
% 
%             disp(['total vol'])
%             approximate_volume, obj.d_region_volume, approximate_volume/obj.d_region_volume

        end
        
        % ======================================================================
        %> @brief Plot true geometry outline.
        % ======================================================================
        function plot_pin(obj)
            P = obj.d_pitch;
            % bounding box
            Xa = [ 0.0;   P;   P; 0.0; 0.0 ];
            Ya = [ 0.0; 0.0;   P;   P; 0.0 ];
            plot (Xa ,Ya ,'k','LineWidth',3);
            hold on
            for i = obj.d_number_radii:-1:1
                t = 0:(2*pi/5000):(2*pi);
                R = obj.d_radii(i);
                plot (P/2+R*cos(t), P/2+R*sin(t), 'k', 'LineWidth', 2);      
            end
            axis([0 P 0 P])
            axis square
            hold off   
        end
        
        % ======================================================================
        %> @brief Plot the mesh.
        % ======================================================================
        function plot_mesh(obj)
            % plot the region map if it exists
            plot_mesh_map(obj, 'REGION')
            hold on
            % Set the mesh map to be semitransparent
            alpha(0.4)
            % Plot true geometry
            plot_pin(obj)
            hold off
        end
        
        % ======================================================================
        %> @brief Plot the tracks (180 degrees).
        % ======================================================================
        function plot_tracks(obj, full)
            if nargin == 1
                full = 0;
            end
            P = obj.d_pitch;
            % bounding box
            Xa = [ 0.0;   P;   P; 0.0; 0.0 ];
            Ya = [ 0.0; 0.0;   P;   P; 0.0 ];
            plot (Xa ,Ya ,'k','LineWidth',3);
            hold on
            % Plot true geometry
            plot_pin(obj)
            hold on  
            q = obj.d_quadrature;
            % Plot the tracks
            for m = 1:number_angles_octant(q)
                hold on
                I = obj.d_enter{m};
                F = obj.d_exit{m};
                col = [rand rand rand];
                for k = 1:number_tracks(q, m)
                    plot ([I(k ,1) ; F(k ,1) ],...
                          [I(k ,2) ; F(k ,2) ],'LineWidth',2,'Color',col);
                end
                if full
                I(:, 1) = pitch(obj)-I(:, 1);
                F(:, 1) = pitch(obj)-F(:, 1);
                col = [rand rand rand];
                for k = 1:number_tracks(q, m)
                    plot ([I(k ,1) ; F(k ,1) ],...
                          [I(k ,2) ; F(k ,2) ],'LineWidth',2,'Color',col);
                end
                end
                axis ([0 pitch(obj) 0 pitch(obj)])
                axis square ;
                grid on;
                hold off;
            end
        end
        
        % ======================================================================
        %> @brief Verify the tracks.
        % ======================================================================
        function verify_tracks(obj)
            % Plot true geometry
            %plot_tracks(obj)
            alpha(0.4)
            hold on  
            q = obj.d_quadrature;
            for m = 1:number_angles_octant(q)
                cos_phi = cos(phi(q, 1, m));
                sin_phi = sin(phi(q, 1, m));
                for t = 1:number_tracks(q, m)
                    x_i = obj.d_enter{m}(t, 1); % incident
                    y_i = obj.d_enter{m}(t, 2); %   x and y
                    for i = 1:obj.d_number_segments(m, t)
                        x_o = x_i + cos_phi * obj.d_segment_length{m}{t}(i);
                        y_o = y_i + sin_phi * obj.d_segment_length{m}{t}(i);
                        reg = obj.d_segment_region{m}{t}(i);
                        plot([x_i x_o],[y_i y_o],'--','LineWidth', ...
                            2 ,'Color',col(obj, reg));
                        x_i = x_o; 
                        y_i = y_o;
                    end
                end
            end
            axis ([0 pitch(obj) 0 pitch(obj)])
            axis square ;
            grid on;
            hold off;
        end
        
        
        function p = pitch(obj)
            p = obj.d_pitch;
        end
        
        function r = radii(obj)
            r = obj.d_radii; 
        end
        
        function bool = tracked(obj)
            bool = obj.d_tracked;
        end
        
        function m = region_mat_map(obj)
            m = obj.d_region_mat_map;
        end

        function v = region_volume(obj)
            v = obj.d_region_volume;
        end   
        
        function n = number_regions(obj)
            n = obj.d_number_regions; 
        end

        function n = number_segments(obj, m, t)
            n = obj.d_number_segments(m, t);
        end
        
        function l = segment_length(obj, m, t, i)
            l = obj.d_segment_length{m}{t}(i);
        end        
        
        function r = segment_region(obj, m, t, i)
            r = obj.d_segment_region{m}{t}(i); 
        end    
        
        function s = space(obj, m)
            s = obj.d_space(m);
        end
        
    end
    
    methods (Access = private)
        
        % ======================================================================
        %> @brief Determine in what region a mesh center resides.
        %
        %>
        %> @param  i 	Horizontal index.
        %> @param  j 	Vertical index.
        %> @return      Region index.
        % ======================================================================
        function r = find_region(obj, i, j)
            x = (i - 0.5) * obj.d_dx(1);
            y = (j - 0.5) * obj.d_dx(1);
            % Loop through the radii.  If I'm in there, that's where I live.
            r = obj.d_number_regions;
            hp = 0.5*obj.d_pitch;
            for p = 1:obj.d_number_radii
               if sqrt((x-hp)^2 + (y-hp)^2) < obj.d_radii(p)
                   r = p;
                   break;
               end
            end     
        end
        
        function color = col(obj, g)
            % this function is a hard-coded color map for the different
            % flux groups.  I've accounted for up to 8 groups.
            % set(ur,'Color',[1 0.7 0.2],'LineWidth',2);
            switch g
                case 1
                    color = [0.0 0.0 1.0]; % blue
                case 2
                    color = [0.0 0.8 0.2]; % nice green
                case 3
                    color = [1.0 0.0 0.0]; % red
                case 4
                    color = [0.4 0.0 0.6]; % purple
                case 5
                    color = [0.9 0.4 0.0]; % orange
                case 6
                    color = [0.5 0.2 0.0]; % brown
                case 7
                    color = [0.0 0.8 0.6]; % turquoise
                case 8
                    color = [0.7 0.6 0.0]; % gold
                otherwise
                    color = [rand rand rand]; % random
            end
        end
        
    end
    
    
end