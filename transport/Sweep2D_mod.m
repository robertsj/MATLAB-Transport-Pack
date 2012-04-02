%> @file  Sweep2D_mod.m
%> @brief Sweep2D_mod class definition.
% ==============================================================================
%> @brief Sweep over a 2D mesh.
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
classdef Sweep2D_mod < SweepBase

    properties        
        %> Place holder horizontal edge flux [num_x]
        d_psi_horizontal      
        %
        d_kernel = 1;
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
        %> @return              Instance of the Sweep2D class.
        % ======================================================================
        function obj = Sweep2D_mod(input, mesh, mat, quadrature, boundary, ...
                                   equation)
            
            % Call base class.
            obj = obj@SweepBase(input, mesh, mat, quadrature, boundary, ...
                                equation);

            
            % Initialize place holders.  Currently, we sweep first along x and
            % then y.  Hence, we need to store a row of horizontal edge fluxes.
            obj.d_psi_horizontal = zeros(number_cells_x(mesh), 1);       
            
        end
        
        % ======================================================================
        %> @brief Sweep the mesh for all angles.
        %
        %> This performs the action of \f$ \mathbf{T}^{-1} \f$ on a given
        %> discrete right hand side vector.  Currently, a single vector
        %> applicable to all angles is provided, since we work only with 
        %> isotropic sources and scattering.
        %>
        %> @param source    Discrete sweep source.
        %> @param g         Group of this problem.
        %>
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
            for o = 1:4
                
                % Get incident boundary fluxes.  These will be written.
                psi_v = get_psi_v_octant(obj.d_boundary, o, Boundary.IN);
                psi_h = get_psi_h_octant(obj.d_boundary, o, Boundary.IN);
                
                % Setup the equations for this octant.
                setup_octant(obj.d_equation, o);
                
                % Get the constants.
                con_x = get_con_x(obj.d_equation, o);   
                con_y = get_con_y(obj.d_equation, o);   
                
                % Get other coefficients.
                beta = obj.d_equation.d_beta;    
                
                % Get mesh bounds.
                xb = x_bounds(obj, o);
                yb = y_bounds(obj, o);
                
                % Use the MEX kernel to sweep.  3x(+?) times faster.
                if obj.d_kernel == 1
                [phi, psi_v, psi_h] = ...
                    sweep2Dkernel_opt2_mex(phi, psi_v, psi_h, obj.d_nx, yb, xb, sig, ...
                                      con_x, con_y, s, wt, beta);
                else
                % Sweep over all cells.
                for j = yb(1):yb(3):yb(2)
                    
                    psi_v_temp = psi_v(j, :); % psi_v(cells, angles)
                    
                    for i = xb(1):xb(3):xb(2)
                        
                        % Set incident angular flux.
                        %psi_in = [psi_h(i, :); psi_v_temp];
                        
                        % Cardinal index.
                        k = i + (j  -1)*obj.d_nx;
                        
                        coef = 1.0 ./ (sig(i, j) + con_x(i, :) + con_y(j, :));
                        psi_center = coef .* ...
                            (s(k) + con_x(i, :) .* psi_v_temp + ...
                                    con_y(j, :) .* psi_h(i, :) );
                                
                        % Outgoing fluxes        
                        psi_h(i, :) = beta(1)*psi_center + beta(2)*psi_h(i, :);
                        psi_v_temp  = beta(1)*psi_center + beta(2)*psi_v_temp;
                        
                        % Inner product of weights with psi.
                        phi(k) = phi(k) + psi_center * wt;
                        
                        % *** Here is where coarse mesh boundary
                        %     information could be dealt with.
                        
                    end
                    
                    % Save outgoing vertical
                    psi_v(j, :) = psi_v_temp;
                    
                end
                end
                % Save outgoing flux vectors.  Note that psi_h has been
                % computed for the last row of the mesh and so is the
                % outgoing boundary flux.
                set_psi_v_octant(obj.d_boundary, o, psi_v, Boundary.OUT);
                set_psi_h_octant(obj.d_boundary, o, psi_h, Boundary.OUT);
                
                
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
        
        function by = y_bounds(obj, o)
            if obj.d_quadrature.octant(o, 2) == 1 
                by = [1 number_cells_y(obj.d_mesh) 1];
            else
                by = [number_cells_y(obj.d_mesh) 1 -1];
            end
        end
        
    end
    
    
    
end
