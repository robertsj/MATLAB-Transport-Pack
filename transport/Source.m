%> @file  Source.m
%> @brief Source class definition.
% ==============================================================================
%> @brief External isotropic volume source.
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
        %> @param   mesh            Mesh.
        %> @param   number_groups   Number of energy groups.
        %> @return                  Instance of the Source class.
        % ======================================================================
        function this = Source(mesh, number_groups)
            this.d_mesh = mesh;
            this.d_number_groups = number_groups;
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
        function this = set_sources(this, sources, source_map)

            DBC.Require('length(sources(:, 1)) == this.d_number_groups');

            % Scale strengths to isotropic using appropriate angular norm.
            this.d_sources = sources * Quadrature.angular_norm(this.d_mesh.DIM);
            
            % ==================================================================
            % MESH
            if meshed(this.d_mesh)
                
                % Add the source map.
                add_mesh_map(this.d_mesh, source_map, 'EXTERNALSOURCE');
                this.d_source_map = mesh_map(this.d_mesh, 'EXTERNALSOURCE');

                % Reshape the source map to a vector.
                src = reshape(this.d_source_map, number_cells(this.d_mesh), 1);

                this.d_source_vec = zeros(number_cells(this.d_mesh),...
                                         this.d_number_groups);

                for cell = 1:number_cells(this.d_mesh)
                    for group = 1:this.d_number_groups
                        this.d_source_vec(cell, group) = ...
                            this.d_sources(group, src(cell));            
                    end
                end
                
            % ==================================================================
            % MOC
            elseif tracked(this.d_mesh)
                
                this.d_source_map = source_map;
                this.d_source_vec = this.d_sources;
                
            else
                error('Invalid spatial mesh! Neither meshed nor tracked.')
            end
            
            this.d_initialized = 1;
            
        end
        
        % ======================================================================
        %> @brief Set external sources on a fine mesh.
        %>
        %> Note, this is a source in terms of *moments*.  The client is 
        %> responsible for applying the moments-to-discrete operator.
        %>
        %> @param sources       Unique sources.
        %> @param source_map 	Fine mesh map of sources.  
        % ======================================================================
        function this = set_sources_mesh(this, sources, source_map)
            this.d_sources    = sources;
            this.d_source_map = source_map;
            this.d_source_vec = zeros(number_cells(this.d_mesh),...
                                     this.d_number_groups);
            for cell = 1:number_cells(this.d_mesh)
                for group = 1:this.d_number_groups
                    this.d_source_vec(cell, group) = ...
                        this.d_sources(group, source_map(cell));
                end
            end
            this.d_initialized = 1;
        end
        
        
        % Getters.
        
        function q = cell_source(this, i, j, g)
            DBC.Require('g > 0 && g <= this.d_number_groups');
            q = this.d_sources(g, this.d_source_map(i, j));
        end
        
        function q = region_source(this, r, g)
            q = this.d_sources(g, this.d_source_map(r));
        end
        
        function q = source(this, g)
            q = this.d_source_vec(:, g);
        end
        
        function y = initialized(this)
            y = this.d_initialized;
        end
        
    end
    
end