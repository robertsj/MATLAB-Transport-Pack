classdef Source < handle
    % External isotropic volume source.
   
    properties
        d_mesh          % Cartesian mesh.
        d_sources       % Unique source spectra (shape + strength).
        d_source_map    % Source fine mesh placement.      
        d_number_groups % Number of groups.
        d_source_vec    % Vector of group-wise isotropic sources
        d_initialized = 0;  % Am I ready?
    end
    
    methods
       
        % ----------------------------------------------------------------------
        % Constructor
        % ----------------------------------------------------------------------
        
        function obj = Source(mesh, number_groups)
            obj.d_mesh = mesh;
            obj.d_number_groups = number_groups;
        end
        
        % ----------------------------------------------------------------------
        % Public Interface
        % ----------------------------------------------------------------------
        
        % Setters.
        
        function obj = set_sources(obj, sources, source_map)
        %  function obj = set_sources(obj, sources, source_map)
            DBC.Require('length(sources(:, 1)) == obj.d_number_groups');
            % DBC.Require('length(material_map(1, :))==length(xfm)')
            % DBC.Require('length(material_map(:, 1))==length(yfm)')

            % Scale strengths to isotropic using appropriate angular norm.
            obj.d_sources = sources * Quadrature.angular_norm(obj.d_mesh.DIM);
            
            % Add the source map.
            obj.d_mesh = add_mesh_map(obj.d_mesh, source_map, 'EXTERNALSOURCE');
            obj.d_source_map = mesh_map(obj.d_mesh, 'EXTERNALSOURCE');
            
            % Reshape the source map to a vector.
            src = reshape(obj.d_source_map, number_cells(obj.d_mesh), 1);
            
            obj.d_source_vec = zeros(number_cells(obj.d_mesh),...
                                     obj.d_number_groups);

            for cell = 1:number_cells(obj.d_mesh)
                for group = 1:obj.d_number_groups
                    obj.d_source_vec(cell, group) = ...
                        obj.d_sources(group, src(cell));            
                end
            end
            
            obj.d_initialized = 1;
        end
        
        % Getters.
        
        function q = cell_source(obj, i, j, g)
            DBC.Require('g > 0 && g <= obj.d_number_groups');
            q = obj.d_sources(g, obj.d_source_map(i, j));
        end
        
        function q = source(obj, g)
            q = obj.d_source_vec(:, g);
        end
        
        function y = initialized(obj)
            y = obj.d_initialized;
        end
        
    end
    
end