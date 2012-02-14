%> @file  SD1D.m
%> @brief SD1D class definition.
% ==============================================================================
%> @brief Diamond difference approximation in one dimension.
%
%> The 1-D discretized transport equation in 
%> @f[
%>     |\mu_n| (\psi_{i,n,out} - \psi_{i,n,in}) + \Delta \Sigma \psi_{i,n}
%>        = \Delta Q_{i, n} \, .
%> @f]
%> To solve this, the cell center flux is related to the incoming and
%> outgoing edge fluxes by a general relation
%> @f[
%>      \psi_{i,n} = \frac{1+\alpha_{i,n}}{2}\psi_{i,n,out} + 
%>                   \frac{1-\alpha_{i,n}}{2}\psi_{i,n,in} \, .
%> @f]
%> Solving for the outgoing flux and inserting in the transport equation
%> yields
%> @f[
%>      \psi_{i,n} = 
%>      \frac{Q_{i,n} + \psi_{i,n,in} \frac{2\mu}{\Delta(\alpha_{i,n}+1)}}
%>           {\frac{2\mu}{\Delta(\alpha_{i,n}+1)} + \Sigma} \, ,
%> @f] 
%> and
%> @f[
%>      \psi_{i,n,out} = \frac{2}{\Delta(\alpha_{i,n}+1)}\psi_{i,n} + 
%>             \frac{\alpha_{i,n}-1}{\alpha_{i,n}+1} \psi_{i,n,in} \, .
%> @f]
%> For the step difference approximation, @f$\alpha = 1 @f$ over
%> all indices.  
%> 
%> In general, we define
%> @f[
%>      con_{x,i,n} = \frac{2\mu}{\Delta(\alpha_{i,n}+1} ,
%> @f]
%> @f[
%>      \beta{1} = \frac{2}{\alpha_{i,n}+1} \, ,
%> @f]
%> and
%> @f[
%>      \beta{2} = \frac{\alpha_{i,n}-1}{\alpha_{i,n}+1} \, .
%> @f]
%> which applies to other discretizations, thus yielding a consistent
%> interface.
%>
%> @sa DD1D
% ==============================================================================
classdef SD1D < Equation
    
    properties
        d_con_x
        d_beta
        d_sig
        d_con_out
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
        %> @return Instance of the SD1D class.
        % ======================================================================
        function obj = SD1D(mesh, mat, quadrature)
            % Call the base class.
            obj = obj@Equation(mesh, mat, quadrature); 
            
            % Presize coefficient vectors.
            obj.d_sig   = zeros(number_cells_x(mesh), 1);
            
            % Get mu's.  By construction, we only need the positive values.
            mu = angle_octant(obj.d_quadrature, 1);

            % Get weight. 
            wt = weight_octant(obj.d_quadrature);
            
            % Get the widths from mesh.
            w = widths(obj.d_mesh);
            
            % Build the angle-dependent constants.
            obj.d_con_x = zeros(length(w{1}), length(wt));
            for i = 1:length(mu)
                obj.d_con_x(:, i) = mu(i)./w{1};
            end
            
            obj.d_beta = [1; 0]; 
            
        end
        
        function c = get_con_x(obj, octant)
            c = obj.d_con_x;
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
            for i = 1:number_cells_x(obj.d_mesh)
                obj.d_sig(i) = sigma_t(obj.d_mat, obj.d_mat_map(i), group);
            end
        end
        
        % ======================================================================
        %> @brief Setup the equations for an octant.
        %>
        %> @param octant    Current octant.
        % ======================================================================
        function obj = setup_octant(obj, octant)
            % DD needs nothing more; everything is done at construction.
        end
        
        % ======================================================================
        %> @brief Setup the equations for an angle.
        %>
        %> @param mu    Cosine with respect to x axis.
        %> @param eta   Cosine with respect to y axis.
        % ======================================================================
        function obj = setup_angle(obj, mu)
            % DD needs nothing more; everything is done at construction.                   
        end
        
        % ======================================================================
        %> @brief Solve for the cell-center and outgoing edge fluxes.
        %>
        %> @param g         Group index
        %> @param psi_in    Incident flux vector
        %> @param s         Cell source
        %> @param i         Cell x index
        %> @return Cell center angular flux and outgoing edge fluxes.
        % ======================================================================
        function [psi_out, psi_center] = solve(obj, g, psi_in, s, i)
% %                     coef       = 1.0 ./ (sig(i) + con_x(i, :));
% %                     psi_center = coef .* (source(i) + con_x(i, :) .* psi_in(1, :));
% %                     psi_out(1, :) = 2*psi_center - psi_in(1, :);            
%             coef = 1.0 ./ (obj.d_sig(i) + obj.d_con_x(i, :));
%             psi_center = coef .* (s + obj.d_con_x(i, :) .* psi_in(1, :));
%             psi_out(1, :) = 2*psi_center - psi_in(1, :);
        end
        
    end
    
end
