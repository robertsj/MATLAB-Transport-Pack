%> @file  Mesh.m
%> @brief Mesh class definition.
% ==============================================================================
%> @brief 2-D Cartesian mesh.
%
%> Finish me.
% ==============================================================================
classdef Mesh < handle
    % Cartesian mesh class.
    
    properties (Constant)
        LEFT    = 1; % left of box
        RIGHT   = 2; % ...
        BOTTOM  = 3;
        TOP     = 4;
        SOUTH   = 5;
        NORTH   = 6;
    end
    
    properties (Constant)
        DIM     = 2;
    end
    
    properties (Access = private)    
        %> x coarse mesh boundaries
        d_xcm   
        %> y coarse mesh boundaries
        d_ycm   
        %> x fine meshes in each x coarse mesh
        d_xfm   
        %> y fine meshes in each y coarse mesh
        d_yfm   
        %> x widths
        d_dx
        %> y widths
        d_dy
        %> Total number of cells
        d_number_cells   
        %> Number of cells in x direction
        d_number_cells_x 
        %> Number of cells in y direction
        d_number_cells_y 
        %> Map container containing a key describing a mesh property and a fine
        %> mesh map defining the property in each cell.  These properties 
        %> include materials, coarse mesh regions (pins, assembly, fuel,
        %> moderator, etc.), and anything else the user wants to edit.
        d_mesh_map       
    end
    
    methods
        
        % ======================================================================
        %> @brief Class constructor
        %>
        %> The mesh is defined in terms of .
        %>
        %> @param xfm           Number of fine mesh per x coarse mesh.
        %> @param yfm           Number of fine mesh per y coarse mesh.
        %> @param xcm           Coarse mesh edges along x axis.
        %> @param ycm           Coarse mesh edges along y axis.
        %> @param material_map	Coarse mesh material map.
        %>
        %> @return Instance of the Mesh class.
        % ======================================================================
        function obj = Mesh(xfm, yfm, xcm, ycm, material_map)
            DBC.Require('length(xfm) > 0') 
            DBC.Require('length(yfm) > 0')
            DBC.Require('length(xcm) > 1')
            DBC.Require('length(ycm) > 1')
            DBC.Require('length(material_map(1, :))==length(xfm)')
            DBC.Require('length(material_map(:, 1))==length(yfm)')
             
            
            obj.d_xfm = xfm;
            obj.d_yfm = yfm;
            obj.d_xcm = xcm;
            obj.d_ycm = ycm;
            obj.d_dx  = zeros(sum(obj.d_xfm), 1);
            obj.d_dy  = zeros(sum(obj.d_yfm), 1);
            obj.d_number_cells_x = sum(obj.d_xfm);
            obj.d_number_cells_y = sum(obj.d_yfm);
            obj.d_number_cells = obj.d_number_cells_x * obj.d_number_cells_y;
            
            % Setup temporary fine mesh material map.
            tmp_material_map = zeros(obj.d_number_cells_x, ...
                                     obj.d_number_cells_y);
            % Flip the material_map.  For user ease, all maps are defined
            % in a natural orientation.  That is, the top row of a matrix
            % corresponds to the physical top row physically, first column
            % mean leftmost physically, etc., all from an overhead point of
            % view.  Computationally, we stick with i->x, j->y.
            material_map = flipud(material_map)';
            
            k  = 0; % temporary place 
            kk = 0; % holders
            
            % Discretize
            for i = 1:length(obj.d_xfm)
                x_range = (k+1):(k+obj.d_xfm(i)); 
                
                % x widths
                obj.d_dx(x_range) = ...
                    (obj.d_xcm(i+1) - obj.d_xcm(i)) / obj.d_xfm(i);
                
                for j = 1:length(obj.d_yfm)
                    y_range = (kk+1):(kk+obj.d_yfm(j)); 
                    
                    % y widths
                    obj.d_dy(y_range) = ...
                        (obj.d_ycm(j+1) - obj.d_ycm(j)) / obj.d_yfm(j);
                    
                    % Assign material region maps
                    tmp_material_map(x_range, y_range) = ...
                        material_map(i, j);
    
                    kk = sum(obj.d_yfm(1:j));
                end
                
                kk = 0;
                k = sum(obj.d_xfm(1:i));
                
            end 
            
            % Insert the material map to the mesh map container.
            obj.d_mesh_map = containers.Map('MATERIAL', ...
                                            tmp_material_map);
        end
       
        
        % ----------------------------------------------------------------------
        % Public Interface
        % ----------------------------------------------------------------------
        
        function obj = add_mesh_map(obj, mesh_map, map_key)
        % function obj = add_mesh_map(obj, mesh_map, map_key)
        %   Add a new mesh map.
        %   
        %   Inputs:
        %     mesh_map -- Map values defined on the *coarse mesh*.  
        %     map_id   -- The integer key given to this map.  These
        % 
            if (isKey(obj.d_mesh_map, map_key))
                disp('Warning: mesh map key is already in use!')
            end
            % Flip the mesh_map
            mesh_map = flipud(mesh_map)';
%             DBC.Require('length(mesh_map(1, :))==length(obj.d_number_cells_x)');
%             DBC.Require('length(mesh_map(1, :))==length(obj.d_number_cells_y)');
            
            tmp_mesh_map = zeros(obj.d_number_cells_x, ...
                                 obj.d_number_cells_y) ;
            
            k  = 0; % temporary place 
            kk = 0; % holders                             
            for i = 1:length(obj.d_xfm)
                x_range = (k+1):(k+obj.d_xfm(i)); 

                for j = 1:length(obj.d_yfm)
                    y_range = (kk+1):(kk+obj.d_yfm(j)); 
                    
                    % Assign coarse mesh vale to the fine mesh
                    tmp_mesh_map(x_range, y_range) = mesh_map(i, j);
    
                    kk = sum(obj.d_yfm(1:j));
                end
                
                kk = 0;
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
        
        function v = plot_mesh_map(obj, mapkey)
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
        
        function v = plot_flux(obj, f)
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