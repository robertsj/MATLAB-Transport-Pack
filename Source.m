%> @file  Source.m
%> @brief Source class definition.
% ==============================================================================
%> @brief External isotropic volume source.
%
%> More here...
% ==============================================================================
classdef Source < handle

    properties
        d_mesh          % Cartesian mesh.
        d_sources       % Unique source spectra (shape + strength).
        d_source_map    % Source fine mesh placement.      
        d_number_groups % Number of groups.
        d_source_vec    % Vector of group-wise isotropic sources
        d_initialized = 0;  % Am I ready?
    end
    
    methods
        
        % ======================================================================
        %> @brief Class constructor
        %>
        %> All boundary flux arrays are sized, and the boundary conditions are
        %> constructed for each surface.
        %>
        %> @param mesh              Mesh.
        %> @param number_groups     Number of energy groups.
        %> @return Instance of the Source class.
        % ======================================================================
        function obj = Source(mesh, number_groups)
            obj.d_mesh = mesh;
            obj.d_number_groups = number_groups;
        end
        
        % ----------------------------------------------------------------------
        % Public Interface
        % ----------------------------------------------------------------------
        
        % ======================================================================
        %> @brief Set external sources.
        %
        %> The sources are defined in columns.  Each column must be equal
        %> to the number of groups.  The map is defined per coarse mesh or
        %> region.  Note, if a region has no source, a zero strength must
        %> be defined.  For the case of MOC, the the source map is a single
        %> vector arranged in the order of the regions.
        %
        %> @param sources    	Unique source spectrum * strength
        %> @param source_map 	Coarse mesh or region map of sources.  
        % ======================================================================
        function obj = set_sources(obj, sources, source_map)

            DBC.Require('length(sources(:, 1)) == obj.d_number_groups');
            % DBC.Require('length(material_map(1, :))==length(xfm)')
            % DBC.Require('length(material_map(:, 1))==length(yfm)')

            % Scale strengths to isotropic using appropriate angular norm.
            obj.d_sources = sources * Quadrature.angular_norm(obj.d_mesh.DIM);
            
            % ==================================================================
            % MESH
            if meshed(obj.d_mesh)
                
                % Add the source map.
                add_mesh_map(obj.d_mesh, source_map, 'EXTERNALSOURCE');
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
                
            % ==================================================================
            % MOC
            elseif tracked(obj.d_mesh)
                
                obj.d_source_map = source_map;
                obj.d_source_vec = obj.d_sources;
                
            else
                error('Invalid spatial mesh! Neither meshed nor tracked.')
            end
            
            obj.d_initialized = 1;
            
        end
        
        
        
        % Getters.
        
        function q = cell_source(obj, i, j, g)
            DBC.Require('g > 0 && g <= obj.d_number_groups');
            q = obj.d_sources(g, obj.d_source_map(i, j));
        end
        
        function q = region_source(obj, r, g)
            q = obj.d_sources(g, obj.d_source_map(r));
        end
        
        function q = source(obj, g)
            q = obj.d_source_vec(:, g);
        end
        
        function y = initialized(obj)
            y = obj.d_initialized;
        end
        
    end
    
end