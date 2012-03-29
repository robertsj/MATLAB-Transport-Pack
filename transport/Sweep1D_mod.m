%> @file  Sweep1D.m
%> @brief Sweep1D class definition.
% ==============================================================================
%> @brief Sweep over a 1D mesh.
%
%> The within-group transport equation is 
%> \f[
%>      \mathbf{T}\psi = Q \, ,
%> \f]
%> where \f$ \mathbf{T} \f$ is the streaming and collision operator and 
%> \f$ Q \f$ is a discrete representation of all source contributions.
%>
%> To invert the operator \f$ \mathbf{T} \f$, we "sweep" over the mesh for all 
%> angles,
%> which gives us updated angular fluxes in each cell.  However, we don't store
%> the angular flux, but rather add its contribution to the scalar flux 
%> directly since storing the angular flux is too expensive.
% ==============================================================================
classdef Sweep1D_mod < SweepBase

    properties (Access = private)
        
    end
    
    methods
        
        % ======================================================================
        %> @brief Class constructor
        %
        %> @param   input       Input database.
        %> @param   mesh      	Cartesian mesh.
        %> @param   mat       	Material definitions.
        %> @param   quadrature  Angular mesh.
        %> @param   boundary    Boundary flux container.
        %> @param   equation    Discretization
        %>
        %> @return              Instance of the Sweep1D class.
        % ======================================================================
        function obj = Sweep1D_mod(input, mesh, mat, quadrature, boundary, ...
                                   equation)
            % Call base class.
            obj = obj@SweepBase(input, mesh, mat, quadrature, boundary, ...
                                equation);
            
        end
        
        % ======================================================================
        %> @brief Sweep the mesh for all angles.
        %
        %> This performs the action of \f$ \mathbf{T}^{-1} \f$ on a given
        %> discrete right hand side vector.  Currently, a single vector
        %> applicable to all angles is provided, since we work only with 
        %> isotropic sources and scattering.
        %>
        %> @param s         Discrete sweep source.
        %> @param g         Group of this problem.
        %> @return          Updated group flux.
        % ======================================================================
        function phi = sweep(obj, s, g)
            
            % Initialize flux update.
            phi = zeros(size(s));
            
            % Get the precomputed cell cross section
            sig = obj.d_equation.d_sig;
            
            % Get weight. 
            wt = weight_octant(obj.d_quadrature);
            
            % Sweep over all octants.
            for o = 1:2   
                
                % Get incident boundary fluxes.  These will be written.
                % psi_v(space, angles)
                psi_o(1, :) = get_psi_v_octant(obj.d_boundary, o, Boundary.IN);
                
                % Setup equation for octant
                setup_octant(obj.d_equation, o);
                
                % Get the constants.
                con_x = get_con_x(obj.d_equation, o);   
                
                % Get other coefficients.
                beta = obj.d_equation.d_beta;                
                
                % Get the sweep bounds.
                xb = x_bounds(obj, o);
                
                % Sweep over all cells.  Angles in this octants are done
                % simultaneosly via vectorization.
                for i = xb(1):xb(3):xb(2)
                    
                    % Set incident angular flux. psi_v has all angles.
                    psi_i(1, :) = psi_o(1, :); 
                    
                    % Solve
                    coef        = 1.0 ./ (sig(i) + con_x(i, :));
                    psi_center  = coef .* (s(i) + con_x(i, :) .* psi_i(1, :));
                    psi_o(1, :) = beta(1)*psi_center + beta(2)*psi_i(1, :);
                    
                    % Inner product of weights with psi
                    phi(i) = phi(i) +  psi_center * wt; 
                    
                end
                
                % Save outgoing flux vectors.
                set_psi_v_octant(obj.d_boundary, o, psi_o, Boundary.OUT);
                
            end

            obj.d_number_sweeps = obj.d_number_sweeps + 1;
            
        end
        
    end
    
    % Implementation
    methods (Access = private)
        
        function bx = x_bounds(obj, o)
        % function bx = x_bounds(obj, octant)
        %   Return starting and ending x indices and the step
            if obj.d_quadrature.octant(o, 1) == 1 
                bx = [1 number_cells_x(obj.d_mesh) 1];
            else
                bx = [number_cells_x(obj.d_mesh) 1  -1];
            end
        end
        
    end
    
    
    
end