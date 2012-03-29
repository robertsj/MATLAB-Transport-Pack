%> @file  Assembly.m
%> @brief Assembly class definition.
% ==============================================================================
%> @brief  Represents reactor assembly of \ref PinCell objects.
%
%> More here.
% ==============================================================================
classdef Assembly < Mesh2D
    
    properties
        
        %> @name Basic Data
        %> @{
        %
        %> Pin cell array
        d_pin
        %> Pin cell map
        d_pin_map
        %> Number of pins in one direction (eg. 17x17)
        d_number_pins_row
        %> Total number of pins (e.g. 289)
        d_number_pins
        %> Assembly pitch
        d_pitch
        %
        %> @}
          
    end
    
    methods
        
        % ======================================================================
        %> @brief Class constructor
        %
        %> @param   pins        Cell array of pin objects.
        %> @param   pin_map     2-D array of pin locations.
        %> @return              Instance of the Assembly class.
        % ======================================================================
        function obj = Assembly(pins, pin_map)
            DBC.Require('length(pins) > 0');
            DBC.Require('max(unique(pin_map)) <= length(pins)');
            obj.d_pin             = pins;
            obj.d_pin_map         = flipud(pin_map)';
            obj.d_number_pins_row = length(obj.d_pin_map(:, 1));
            obj.d_number_pins     = length(length(pins));   
            obj.d_pitch           = obj.d_number_pins_row * pitch(pins{1});
        end
        
        % ======================================================================
        %> @brief Construct a mesh over the assembly.
        %
        %> Using the underlying pin cell meshes, this constructs an assembly
        %> mesh.  The pin cell mesh maps are translated to the
        %> assembly.  More work will be to differentiate pins.
        % ======================================================================
        function obj = meshify(obj)
            
            % Set number of cells.  This *assumes* all pins have the same
            % meshing, though not necessarily uniform.
            xfm = number_cells_x(obj.d_pin{1}); % pin cell fine
            yfm = number_cells_y(obj.d_pin{1}); % mesh counts
            obj.d_number_cells_x = xfm * obj.d_number_pins_row;
            obj.d_number_cells_y = yfm * obj.d_number_pins_row;
            obj.d_number_cells = obj.number_cells_x * obj.number_cells_y;
            obj.d_xfm = ones(obj.d_number_cells_x, 1);
            obj.d_yfm = ones(obj.d_number_cells_x, 1);
            % Compute the widths.
            obj.d_dx = zeros(obj.number_cells_x, 1);
            obj.d_dy = zeros(obj.number_cells_y, 1);
            w = widths(obj.d_pin{1});
            dx = w{1};
            dy = w{2};
            
            tmp_material_map = zeros(obj.number_cells_x, obj.number_cells_y); 
            tmp_region_map   = zeros(obj.number_cells_x, obj.number_cells_y); 
            
            for j = 1:obj.d_number_pins_row
                
                % Fine mesh y range.
                y_range = 1 + yfm * (j - 1) : j * yfm;
                
                % Set fine mesh width in y direction.
                obj.d_dy(y_range) = dy;
                
                for i = 1:obj.d_number_pins_row
                    
                    % Fine mesh x range.
                    x_range = 1 + xfm * (i - 1) : i * xfm;
                    
                    % Set fine mesh width in x direction.
                    if j == 1
                        obj.d_dx(x_range) = dx;
                    end
                    
                    % Cardinal pin index.
                    p = i + (j - 1) * obj.d_number_pins_row;

                    % Get the material and region maps for this pin.
                    pin = obj.d_pin{obj.d_pin_map(i, j)};
                    m = mesh_map(pin, 'MATERIAL');
                    r = mesh_map(pin, 'REGION');
                    
                    % Assign the values.
                    tmp_material_map(x_range, y_range)  = m;
                    tmp_region_map(x_range, y_range)    = r;

                end
                
            end
            
            % Add maps
            obj.d_mesh_map = containers.Map('MATERIAL', tmp_material_map);
            obj.d_mesh_map('REGION') = tmp_region_map;
            
        end % end function meshify
   
        function p = pitch(obj)
            p = obj.d_pitch;
        end
        
        function v = volume(obj)
            v = obj.d_pitch^2;
        end
        
        
    end
    
    methods (Access = private)
       
        
    end
    
end