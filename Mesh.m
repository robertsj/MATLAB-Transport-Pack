%> @file  Mesh.m
%> @brief Mesh class definition.
% ==============================================================================
%> @brief Base Cartesian mesh.
%
%> Finish me.
% ==============================================================================
classdef Mesh < handle
    
    properties (Constant)
        LEFT    = 1; % left of box
        RIGHT   = 2; % ...
        BOTTOM  = 3;
        TOP     = 4;
        SOUTH   = 5;
        NORTH   = 6;
    end
    
    properties (Access = protected)    
        %> x coarse mesh boundaries
        d_xcm = 1;
        %> y coarse mesh boundaries
        d_ycm = 1;
        %> z coarse mesh boundaries
        d_zcm = 1;
        %> x fine meshes in each x coarse mesh
        d_xfm = 1;
        %> y fine meshes in each y coarse mesh
        d_yfm = 1;
        %> z fine meshes in each y coarse mesh
        d_zfm = 1;
        %> x widths
        d_dx = 1;
        %> y widths
        d_dy = 1;
        %> z widths
        d_dz = 1;
        %> Total number of cells
        d_number_cells   = 1;
        %> Number of cells in x direction
        d_number_cells_x = 1;
        %> Number of cells in y direction
        d_number_cells_y = 1;
        %> Number of cells in y direction
        d_number_cells_z = 1;
        %> Map container containing a key describing a mesh property and a fine
        %> mesh map defining the property in each cell.  These properties 
        %> include materials, coarse mesh regions (pins, assembly, fuel,
        %> moderator, etc.), and anything else the user wants to edit.
        d_mesh_map 
        %> Flag indicating I'm meshed.
        d_meshed = 1; 
    end
    
    methods
        
        % ======================================================================
        %> @brief Class constructor
        %> @return Instance of the Mesh class.
        % ======================================================================
        function obj = Mesh()

        end
       
        
        % **********************************************************************
        % Public Interface
        % **********************************************************************
        
        
        % ======================================================================
        %> @brief Get cell array of width vectors.
        %> @return  Cell array of width vectors.
        % ======================================================================
        function w = widths(obj)
            w = {obj.d_dx, obj.d_dy, obj.d_dz};
        end
        
        % ======================================================================
        %> @brief   Width along x-axis for a cell.
        %> @return  Width.
        % ======================================================================
        function v = dx(obj, i)
           v = obj.d_dx(i); 
        end
        
        % ======================================================================
        %> @brief   Width along y-axis for a cell.
        %> @return  Width.
        % ======================================================================
        function v = dy(obj, i)
           v = obj.d_dy(i); 
        end
        
        % ======================================================================
        %> @brief   Width along z-axis for a cell.
        %> @return  Width.
        % ======================================================================
        function v = dz(obj, i)
           v = obj.d_dz(i); 
        end 
        
        % ======================================================================
        %> @brief   Get the total number of cells.
        %> @return  Number of cells.
        % ======================================================================
        function n = number_cells(obj)
            n = obj.d_number_cells;
        end
        
        % ======================================================================
        %> @brief   Get the number of cells along the x axis.
        %> @return  Number of cells.
        % ======================================================================
        function n = number_cells_x(obj)
            n = obj.d_number_cells_x;
        end
        
        % ======================================================================
        %> @brief   Get the number of cells along the y axis.
        %> @return  Number of cells.
        % ======================================================================
        function n = number_cells_y(obj)
            n = obj.d_number_cells_y;
        end
        
        % ======================================================================
        %> @brief   Get the number of cells along the z axis.
        %> @return  Number of cells.
        % ======================================================================
        function n = number_cells_z(obj)
            n = obj.d_number_cells_z;
        end
        
        % ======================================================================
        %> @brief   Get the number of cells along the an axis.
        %> @return  Number of cells.
        % ======================================================================
        function d = dim(obj, i)
            DBC.Require('i > 0 && i < 4');
            if i == 1
                d = obj.d_number_cells_x;
            elseif i == 2
                d = obj.d_number_cells_y;
            else
                d = obj.d_number_cells_z;
            end
        end
        
        
        % ======================================================================
        %> @brief   Returns the cardinal index for i, j, and k
        %
        %> For efficiency, the client may want to hardcode the indexing into
        %> their routine, as this is suboptimal within a loop.
        %>
        %> @return  Index.
        % ======================================================================
        function idx = index(obj, i, j, k)
            if nargin < 4
                k = 1;
            end
            if nargin < 3
                j = 1;
            end   
            idx = i + ...
                  (j - 1)*obj.d_number_cells_x + ...
                  (k - 1)*obj.d_number_cells_x*obj.d_number_cells_y;   
            % Postcondition.
            DBC.Ensure('idx > 0 && idx <= number_cells(obj)')
        end
        
        % ======================================================================
        %> @brief Get a mesh map.
        %
        %> @param   mesh_key	A mesh map key.
        %> @return              The map.
        % ======================================================================
        function m = mesh_map(obj, map_key)
            m = obj.d_mesh_map(map_key);
        end
        
        function bool = meshed(obj)
            bool = obj.d_meshed; 
        end
        
        
        
        % **********************************************************************
        % Abstract Interface
        % **********************************************************************
        
        
        % ======================================================================
        %> @brief Add a map to the mesh.
        %
        %> @param mesh_map      Map on the coarse mesh.
        %> @param map_key       Key (descriptor) for this map.
        % ======================================================================
        obj = add_mesh_map(obj, mesh_map, map_key)
 
        % ======================================================================
        %> @brief Plot a mesh map on the mesh.
        %> @param map_key 	A mesh map key.
        % ======================================================================
        plot_mesh_map(obj, map_key)

        % ======================================================================
        %> @brief Plot a given flux or flux-shaped vector on the mesh.
        %> @param f 	A vector of values that live on the mesh.
        % ======================================================================
        plot_flux(obj, f)
        
    end
    
end