%> @file  Sweep2D.m
%> @brief Sweep2D class definition.
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
classdef Sweep2D < handle

    properties
        d_input                 % User input.
        d_mesh                  % Cartesian mesh.
        d_mat                   % Materials.
        d_quadrature            % Angular mesh.
        %
        d_equation              % Spatial discretization
        %
        d_alpha                 % Weighted diamond-difference parameters.
        %
        d_psi_horizontal        % place holder horizontal flux [num_x]
        %
        d_boundary
        %
        d_nx
        d_ny
        
    end
    
    methods
        
        % ======================================================================
        %> @brief Class constructor
        %>
        %> 
        %>
        %> @param input        	Number of fine mesh per x coarse mesh.
        %> @param mesh         	Number of fine mesh per y coarse mesh.
        %> @param mat           Coarse mesh edges along x axis.
        %> @param quadrature   	Coarse mesh edges along y axis.
        %> @param boundary      Coarse mesh material map.
        %>
        %> @return Instance of the Mesh class.
        % ======================================================================
        function obj = Sweep2D(input, mesh, mat, quadrature, boundary, equation)
            
            obj.d_input      = input;
            obj.d_mesh       = mesh;
            obj.d_mat        = mat;
            obj.d_quadrature = quadrature;
            obj.d_boundary   = boundary;
            obj.d_equation   = equation;
            
            % Store some things from the mesh.
            obj.d_nx = number_cells_x(mesh);
            obj.d_ny = number_cells_y(mesh);
            
            % Initialize place holders.  Currently, we sweep first along x and
            % then y.  Hence, we need to store a row of horizontal edge fluxes.
            obj.d_psi_horizontal = zeros(number_cells_x(obj.d_mesh), 1);       
            
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
        function phi = sweep(obj, source, g)
            
            phi = zeros(size(source));
            
            % Sweep over all octants.
            for o = 1:4
                
                %disp([' OCTANT = ', num2str(o)])
                
                % Setup the equations for this octant.
                setup_octant(obj.d_equation, o);
                
                % Sweep over all angles in the octant.
                for a = 1:number_angles_octant(obj.d_quadrature)
                    

                    % Get incident boundary fluxes.  These will be written.
                    psi_v = get_psi_v(obj.d_boundary, o, a);
                    psi_h = get_psi_h(obj.d_boundary, o, a);
                    
                    % Get cosines.
                    [mu, eta] = angle(obj.d_quadrature, o, a);
                    
                    % Get weight. 
                    w = weight(obj.d_quadrature, a);
                    
                    % Get mesh bounds.
                    xb = x_bounds(obj, o);
                    yb = y_bounds(obj, o);
                    
                    % Setup equation for this angle.
                    setup_angle(obj.d_equation, mu, eta);
                    
                    % Sweep over all cells.
                    for j = yb(1):yb(3):yb(2)
                        
                        psi_v_temp = psi_v(j);
                        
                        for i = xb(1):xb(3):xb(2)
                            
                            % Set incident angular flux.
                            psi_in = [psi_h(i); psi_v_temp];
                            
                            % Cardinal index.
                            k = i + (j-1)*obj.d_nx;% index(obj.d_mesh, i, j);
                            
                            % Solve for this cell.
                            [psi_out, psi_center] = ...
                                solve(obj.d_equation, ...
                                    psi_in, source(k), i, j, g);
                            phi(k) = phi(k) + w*psi_center;
                            % Save the outgoing angular flux.
                            psi_h(i) = psi_out(1);
                            psi_v_temp = psi_out(2);
                        
                            % *** Here is where coarse mesh boundary
                            %     information could be dealt with.
                            
                        end
                        
                        % Save outgoing vertical
                        psi_v(j) = psi_v_temp;
                        
                        
                    
                    end
                    % Save outgoing flux vectors.  Note that psi_h has been
                    % computed for the last row of the mesh and so is the
                    % outgoing boundary flux.
                    set_psi_v(obj.d_boundary, o, a, psi_v);
                    set_psi_h(obj.d_boundary, o, a, psi_h);
  
                end
                
            end
       
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