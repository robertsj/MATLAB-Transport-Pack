%> @file  PinCell.m
%> @brief PinCell class definition.
% ==============================================================================
%> @brief  Represents reactor assembly pin in two dimensions.
%
%> More here.
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
        
        % Remember we get the following mesh data:
        %   d_xcm = 1;
        %   d_ycm = 1;
        %   d_zcm = 1;
        %   d_xfm = 1;
        %   d_yfm = 1;
        %   d_zfm = 1;
        %   d_dx = 1;
        %   d_dy = 1;
        %   d_dz = 1;
        %   d_number_cells   = 1;
        %   d_number_cells_x = 1;
        %   d_number_cells_y = 1;
        %   d_number_cells_z = 1;
        %   d_mesh_map       
        
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
        %> @param  pitch 	Pin cell width
        %> @param  radii 	Fuel pin radii, if any.
        %> @param  mat_map	Material identifiers for each region.
        %> @return Instance of the PinCell class.
        % ======================================================================
        function obj = PinCell(pitch, radii, mat_map)
            % Precondition
            DBC.Require('length(radii) + 1 == length(mat_map)');
            
            obj.d_pitch          = pitch;
            obj.d_center         = [0.5*pitch 0.5*pitch];
            obj.d_radii          = radii;
            obj.d_number_radii   = length(radii);
            obj.d_number_regions = length(radii) + 1;
            
            % Don't allow radii >= to the pitch for simplicity.
            if max(radii) >= pitch/2
                error('Radii must be smaller than half pitch')
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
            
        end
        
        function plot_pin(obj)
            
            P = obj.d_pitch;

            % bounding box
            Xa = [ 0.0;   P;   P; 0.0; 0.0 ];
            Ya = [ 0.0; 0.0;   P;   P; 0.0 ];
            plot (Xa ,Ya ,'k','LineWidth',4);
            hold on;
            
            % plot the region map if it exists
            if ~isempty(obj.d_dx)
                plot_mesh_map(obj, 'REGION')
            end
            
            % Set the real geometry to to be a transparent background
            alpha(0.4)
            
            % plotting pins
            for i = obj.d_number_radii:-1:1
                t = 0:(2*pi/1000):(2*pi);
                R = obj.d_radii(i);
                plot (P/2+R*cos(t), P/2+R*sin(t), 'k', 'LineWidth', 3);      
            end
            axis([0 P 0 P])
            axis square
            hold off
        end
        
        function p = pitch(obj)
            p = obj.d_pitch;
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
        
        function color = col(g)
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