%> @file  Core.m
%> @brief Core class definition.
% ==============================================================================
%> @brief  Represents reactor Core.
%
%> More here.
% ==============================================================================
classdef Core < Mesh2D
    
    properties
        %> Cell array of assembly objects.
        d_assembly
        %> Array of assembly locations.
        d_assembly_map
        %> Number of assemblies per row.
        d_number_assemblies_row
        %> Total number of assemblies.
        d_number_assemblies    
    end
    
    methods
        
        % ======================================================================
        %> @brief Class constructor
        %
        %> @param   assemblies      Cell array of assembly objects.
        %> @param   assembly_map	2-D array of assembly locations.
        %> @return                  Instance of the Core class.
        % ======================================================================
        function obj = Core(assemblies, assembly_map)
            DBC.Require('length(assemblies) > 0');
            DBC.Require( ...
                'length(reshape(assembly_map,[],1)) == length(assemblies)');
            obj.d_assembly              = assemblies;
            obj.d_assembly_map          = flipud(assembly_map)';
            obj.d_number_assemblies_row = length(obj.d_assembly_map(:, 1));
            obj.d_number_assemblies     = length(length(assemblies));   
        end
        
        % ======================================================================
        %> @brief Construct a mesh over the core.
        %
        %> Using the underlying assembly meshes, this constructs 
        %> a core
        %> mesh.  The assembly mesh maps are translated to the
        %> assembly.  More work will be to differentiate pins.
        % ======================================================================
        function obj = meshify(obj)
            
            % Set number of cells.  This *assumes* all assemblies have the 
            % same meshing, though not necessarily uniform.
            xfm = number_cells_x(obj.d_assembly{1}); % assembly fine
            yfm = number_cells_y(obj.d_assembly{1}); % mesh counts
            
            % We're also assuming square cores for now.
            obj.d_number_cells_x = xfm * obj.d_number_assemblies_row;
            obj.d_number_cells_y = yfm * obj.d_number_assemblies_row;
            obj.d_number_cells   = obj.number_cells_x * obj.number_cells_y;
            
            % Compute the widths.
            obj.d_dx = zeros(obj.number_cells_x, 1);
            obj.d_dy = zeros(obj.number_cells_y, 1);
            w = widths(obj.d_assembly{1});
            dx = w{1};
            dy = w{2};
            
            tmp_material_map = zeros(obj.number_cells_x, obj.number_cells_y); 
            tmp_region_map   = zeros(obj.number_cells_x, obj.number_cells_y); 
            
            for j = 1:obj.d_number_assemblies_row
                
                % Fine mesh y range.
                y_range = 1 + yfm * (j - 1) : j * yfm;
                
                % Set fine mesh width in y direction.
                obj.d_dy(y_range) = dy;
                
                for i = 1:obj.d_number_assemblies_row
                    
                    % Fine mesh x range.
                    x_range = 1 + xfm * (i - 1) : i * xfm;
                    
                    % Set fine mesh width in x direction.
                    if j == 1
                        obj.d_dx(x_range) = dx;
                    end
                    
                    % Cardinal pin index.
                    p = i + (j - 1) * obj.d_number_assemblies_row;

                    % Get the material and region maps for this pin.
                    pin = obj.d_assembly{obj.d_assembly_map(i, j)};
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
        
    end
    
end