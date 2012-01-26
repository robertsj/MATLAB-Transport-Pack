%> @file  SweepMOC.m
%> @brief SweepMOC class definition.
% ==============================================================================
%> @brief Sweep over a tracked 2D domain.
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
classdef SweepMOC < handle

    properties
        
        %> User input.
        d_input      
        %> Tracking information.
        d_track
        %> Materials.
        d_mat                   
        %> Angular mesh.
        d_quadrature            
        %> Spatial discretization
        d_equation              
        %> Weighted diamond-difference parameters.
        d_alpha                 
        %
        %d_psi_horizontal        % place holder horizontal flux [num_x]
        %> Boundary fluxes
        d_boundary
        %
        d_nx
        %d_ny
        
    end
    
    methods
        
        % ======================================================================
        %> @brief Class constructor
        %
        %> @param input        	User input database.
        %> @param mesh         	Tracked domain.
        %> @param mat           Materials definitions.
        %> @param quadrature   	Angular quadrature.
        %> @param boundary      Boundary fluxes.
        %> @param equation      Spatial discretization
        %> @return              Instance of the SweepMOC class.
        % ======================================================================
        function obj = SweepMOC(input, mesh, mat, quadrature, boundary, equation)
            
            obj.d_input      = input;
            obj.d_track      = mesh;
            obj.d_mat        = mat;
            obj.d_quadrature = quadrature;
            obj.d_boundary   = boundary;
            obj.d_equation   = equation;
            
            % Store some things from the mesh.
            obj.d_nx = number_cells_x(mesh);
            
            % Initialize place holders.  Currently, we sweep first along x and
            % then y.  Hence, we need to store a row of horizontal edge fluxes.
            %obj.d_psi_horizontal = zeros(number_cells_x(obj.d_track), 1);       
            
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
        %> @return          Updated group flux.
        % ======================================================================
        function phi = sweep(obj, source, g)
            
            % 
            IN  = 1;
            OUT = 2;
            
            % Temporary region flux
            phi = zeros(size(source));

            % Material map
            mat_map = region_mat_map(obj.d_track);

            % Inverse of region volumes.
            inv_volume = 1 ./ region_volume(obj.d_track);
            vol = 0*inv_volume;
            % Single track segment incident and exiting flux
            psi_out = zeros(number_polar(obj.d_quadrature), 1);
            psi_in  = zeros(number_polar(obj.d_quadrature), 1);
            
            % Sweep over all octants.
            for o = 1:4

                % Setup the equations for this octant.
                setup_octant(obj.d_equation, o);
                
                % Sweep over all azimuthal angles in the octant.
                for a = 1:number_angles_octant(obj.d_quadrature)   
                    
                    % Get azimuthal weight. 
                    w_a = weight_phi(obj.d_quadrature, a);
%                     w_a = 1;
                    
                    % Get incident boundary fluxes for this angle.  These
                    % are in all polar angles
                    psi = get_psi(obj.d_boundary, o, a, IN);

                    % Track spacing
                    spacing = space(obj.d_track, a);

                    % Sweep along all tracks
                    for t = 1:number_tracks(obj.d_quadrature, a)

                        % Set the incident flux.
                        psi_out(:, 1) = psi(:, t);

                        % Sweep over all segments on this track
                        for i = 1:number_segments(obj.d_track, a, t)

                            % Set incident angular flux.  (polar, track)
                            psi_in(:, 1) = psi_out(:, 1);

                            % Region
                            r = segment_region(obj.d_track, a, t, i);
                            
                            % Segment length
                            len = segment_length(obj.d_track, a, t, i);
                            
                            % Cross section
                            sigma = sigma_t(obj.d_mat, mat_map(r), g);
                            
                            % Solve for this cell over all polar angles.
                            for p = 1:number_polar(obj.d_quadrature)
                                
                                % Scaled segment length
                                tau = len / mu(obj.d_quadrature, p);
                                
                                % Solve for outgoing and average segment flux.
                                [psi_out(p, 1), psi_avg] = ... 
                                    solve(obj.d_equation, ...
                                          psi_in(p, 1), ...
                                          source(r), ...
                                          sigma, ...
                                          tau);

                                % Polar weight.
                                w_p = weight_mu(obj.d_quadrature, p);
%                                 w_p = 1;
                               % psi_out(p, 1) = 1; psi_avg =1;
                                % Integration of scalar flux.  When all
                                % sweeps are done, this is equivalent to
                                % e.g. Eq. 3.347 in Hebert.
                                phi(r) = phi(r) +  w_a * w_p * ...
                                	psi_avg * spacing * len * inv_volume(r);
                                

                                % 4*pi*V(r) for debugging.
                                vol(r) = vol(r) +  psi_avg;% w_a * w_p * spacing * len;
                                % Note, this *assumes* the spacing is
                                % corrected!!

                            end % polar

                        end % segment
                        
                        % Save the outgoing track flux.
                        psi(:, t) = psi_out(:, 1);

                    end % track
                    
                    % Save outgoing flux vectors for this angle.
                    set_psi(obj.d_boundary, o, a, psi, OUT);
                
                end % azimuth
                
     
            end
            disp('summed psi...')
            vol
%             disp('---------')
        end
        
    end
    
    % Implementation
    methods (Access = private)
        
        function bx = x_bounds(obj, o)
        % function bx = x_bounds(obj, octant)
        %   Return starting and ending x indices and the step
            if obj.d_quadrature.octant(o, 1) == 1 
                bx = [1 number_x(obj.d_quadrature) 1];
            else
                bx = [number_x(obj.d_quadrature) 1  -1];
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