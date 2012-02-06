%> @file  DD2D.m
%> @brief DD2D class definition.
% ==============================================================================
%> @brief Diamond difference approximation in two dimensions.
%
%> 
%> Finish me.
% ==============================================================================
classdef DD2D < Equation
    
    properties
        d_con_x
        d_con_y
        d_sig_all
        d_sig
        d_beta
    end
    
    properties (Constant)
%         HORZ = 1;
%         VERT = 2;
    end
    
    methods
        
        % ======================================================================
        %> @brief Class constructor
        %
        %> Set the mesh and material.
        %>
        %> @param mesh        	Problem mesh.
        %> @param mat          	Material definitions.
        %>
        %> @return Instance of the Equation class.
        % ======================================================================
        function obj = DD2D(mesh, mat, quadrature)
            % Call the base class.
            obj = obj@Equation(mesh, mat, quadrature); 
            
            % Presize coefficient vectors.
            obj.d_sig_all   = zeros(number_cells_x(mesh), ...
                number_cells_y(mesh), number_groups(mat));
            
            % Get mu's.  By construction, we only need the positive values.
            mu  = angle_octant(obj.d_quadrature, 1);
            eta = angle_octant(obj.d_quadrature, 2);

            % Get weight. 
            wt = weight_octant(obj.d_quadrature);
            
            % Get the widths from mesh.
            w = widths(obj.d_mesh);
            
            % Build the angle-dependent constants.  (For DD only right now)
            obj.d_con_x = zeros(length(w{1}), length(wt));
            for i = 1:length(mu)
                obj.d_con_x(:, i) = 2*mu(i)./w{1};
            end
            obj.d_con_y = zeros(length(w{1}), length(wt));
            for i = 1:length(eta)
                obj.d_con_y(:, i) = 2*eta(i)./w{1};
            end
            obj.d_beta = [2; -1]; 
            
            %
            for g = 1:number_groups(mat)
                for j = 1:number_cells_y(obj.d_mesh)
                    for i = 1:number_cells_x(obj.d_mesh)
                        obj.d_sig_all(i, j, g) = ...
                            sigma_t(obj.d_mat, obj.d_mat_map(i, j), g);
                    end
                end
            end
            
        end
        
        function c = get_con_x(obj, octant)
            c = obj.d_con_x;
        end    
        function c = get_con_y(obj, octant)
            c = obj.d_con_y;
        end           
        
        % ======================================================================
        %> @brief Setup the equations for a group.
        %
        %> Here, we'll go through the grid and produce a fine mesh matrix
        %> of total cross-sections.  This isn't the best thing for memory,
        %> but it cuts down a lot on the time within solve.
        %>
        %> @param group     Current group.
        % ======================================================================
        function obj = setup_group(obj, group)
            % 
%             for j = 1:number_cells_y(obj.d_mesh)
%                 for i = 1:number_cells_x(obj.d_mesh)
%                     obj.d_sig(i, j) = ...
%                         sigma_t(obj.d_mat, obj.d_mat_map(i, j), group);
%                 end
%             end
            obj.d_sig = obj.d_sig_all(:, :, group);
        end
        
        % ======================================================================
        %> @brief Setup the equations for an octant.
        %>
        %> @param octant    Current octant.
        % ======================================================================
        function obj = setup_octant(obj, octant)
            % Nothing here for now.   
        end
        
        % ======================================================================
        %> @brief Setup the equations for an angle.
        %>
        %> @param mu    Cosine with respect to x axis.
        %> @param eta   Cosine with respect to y axis.
        %> @param xi    Cosine with respect to y axis.
        % ======================================================================
        function obj = setup_angle(obj, mu, eta, xi)
            % Get the widths from mesh.
            w = widths(obj.d_mesh);
            % Build the two constants.
            obj.d_con_x = 2*abs(mu)./w{1};
            obj.d_con_y = 2*abs(eta)./w{2};    
        end
        
        % ======================================================================
        %> @brief Solve for the cell-center and outgoing edge fluxes.
        %
        %> @param g         Group index
        %> @param psi_in    Incident flux vector, [horz=1, vert=2]
        %> @param s         Cell source
        %> @param i         Cell x index
        %> @param j         Cell y index
        %> @param k         Cell z index (NOT PRESENT)
        %>
        %> @return Cell center angular flux and outgoing edge fluxes.
        % ======================================================================
        function [psi_out, psi_center] = solve(obj, g, psi_in, s, i, j, k)
            %sig  = obj.d_sig(i, j);
            coef = 1.0 / (obj.d_sig(i, j) + obj.d_con_x(i) + obj.d_con_y(j));
            psi_center = coef * (s + obj.d_con_x(i) * psi_in(2) + ...
                                     obj.d_con_y(j) * psi_in(1) );
            psi_out(1) = 2*psi_center - psi_in(1);
            psi_out(2) = 2*psi_center - psi_in(2);   
        end
        
    end
    
end