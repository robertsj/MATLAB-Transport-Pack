%> @file  Mesh1D.m
%> @brief Mesh class definition.
% ==============================================================================
%> @brief 1-D Cartesian mesh.
%
%> Finish me.
% ==============================================================================
classdef Mesh1D < Mesh

    properties (Constant)
        DIM = 1;
    end

    methods
        
        % ======================================================================
        %> @brief Class constructor
        %
        %> @param xfm           Number of fine mesh per x coarse mesh.
        %> @param xcm           Coarse mesh edges along x axis.
        %> @param material_map	Coarse mesh material map.
        %>
        %> @return Instance of the Mesh1D class.
        % ======================================================================
        function obj = Mesh1D(xfm, xcm, material_map)
            
            DBC.Require('length(xfm) > 0') 
            DBC.Require('length(xcm) > 1')
            DBC.Require('length(material_map)==length(xfm)')
             
            obj.d_xfm = xfm;
            obj.d_xcm = xcm;
            obj.d_dx  = zeros(sum(obj.d_xfm), 1);
            obj.d_number_cells_x = sum(obj.d_xfm);
            obj.d_number_cells = obj.d_number_cells_x;
            
            % Setup temporary fine mesh material map.
            tmp_material_map = zeros(obj.d_number_cells, 1);

            k  = 0; % temporary place holder

            % Discretize
            for i = 1:length(obj.d_xfm)
                
                % Fine mesh indices for this coarse mesh.
                x_range = (k+1):(k+obj.d_xfm(i)); 
                
                % Compute the fine mesh widths.
                obj.d_dx(x_range) = ...
                    (obj.d_xcm(i+1) - obj.d_xcm(i)) / obj.d_xfm(i);
                
                % Assign material region maps.
                tmp_material_map(x_range, 1) = material_map(i);

                k = sum(obj.d_xfm(1:i));
                
            end 
            
            % Insert the material map to the mesh map container.
            obj.d_mesh_map = containers.Map('MATERIAL', ...
                                            tmp_material_map);
        end
       
        
        % ----------------------------------------------------------------------
        % Public Interface
        % ----------------------------------------------------------------------
        
        
        % ======================================================================
        %> @brief Add a map to the mesh.
        %
        %> @param mesh_map      Map on the coarse mesh.
        %> @param map_key       Key (descriptor) for this map.
        % ======================================================================
        function obj = add_mesh_map(obj, mesh_map, map_key)

            if (isKey(obj.d_mesh_map, map_key))
                disp('Warning: mesh map key is already in use!')
            end

            tmp_mesh_map = zeros(obj.d_number_cells_x, 1);
            
            k  = 0; % temporary place holder
            
            for i = 1:length(obj.d_xfm)
                
                % Fine mesh indices for this coarse mesh.
                x_range = (k+1):(k+obj.d_xfm(i)); 

                % Assign coarse mesh vale to the fine mesh
                tmp_mesh_map(x_range, 1) = mesh_map(i);
 
                k = sum(obj.d_xfm(1:i));
                
            end   
            
            % Insert the mesh map to the mesh map container.
            obj.d_mesh_map(map_key) = tmp_mesh_map;
        end
        
        % Getters
        
        function w = widths(obj)
        % function w = widths(obj, i, j)
        %   Returns a two element array of cell widths
            w = {obj.d_dx, obj.d_dy};
        end
        
        function v = dx(obj, i)
           v = obj.d_dx(i); 
        end
        
        function v = dy(obj, i)
           v = obj.d_dy(i); 
        end
        
        function n = number_cells(obj)
            n = obj.d_number_cells;
        end
        
        function n = number_cells_x(obj)
            n = obj.d_number_cells_x;
        end
        
        function n = number_cells_y(obj)
            n = obj.d_number_cells_y;
        end
        
        function d = dim(obj, i)
            DBC.Require('i > 0 && i < 3');
            if i == 1
                d = obj.d_number_cells_x;
            else
                d = obj.d_number_cells_y;
            end
        end
        
        function k = index(obj, i, j)
        % function y = index(i, j)
        %   Returns the cardinal index for i and j
            k = i + (j - 1)*obj.d_number_cells_x;
        %    DBC.Ensure('k > 0 && k <= number_cells(obj)')
        end
        
        function m = mesh_map(obj, map_key)
            m = obj.d_mesh_map(map_key);
        end
        
        function plot_mesh_map(obj, mapkey)
            [x, y] = mesh_axes(obj);
            [X, Y] = meshgrid(x, y);
            m = mesh_map(obj, mapkey);
            M = zeros(length(x), length(y));
            M(1:end-1, 1:end-1) = m;
            pcolor(X, Y, M)
            xlabel('x [cm]')
            ylabel('y [cm]')
            title(mapkey), colormap('jet'), colorbar
        end
        
        function plot_flux(obj, f)
            [x, y] = mesh_axes(obj);
            [X, Y] = meshgrid(x, y);
            m = reshape(f, obj.d_number_cells_x, obj.d_number_cells_y);
            M = zeros(length(x), length(y));
            M(1:end-1, 1:end-1) = m;
            pcolor(X, Y, M)
            xlabel('x [cm]')
            ylabel('y [cm]')
            colormap('jet'), colorbar
        end
        
        function [x, y] = mesh_axes(obj)
            x = zeros(obj.d_number_cells_x + 1, 1);
            y = zeros(obj.d_number_cells_y + 1, 1);
            for i = 1:obj.d_number_cells_x
                x(i + 1) = x(i) + obj.d_dx(i);
            end
            for i = 1:obj.d_number_cells_y
                y(i + 1) = y(i) + obj.d_dy(i);
            end
        end
        
    end
    
end