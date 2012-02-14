%> @file  Reflective.m
%> @brief Reflective boundary condition class definition.
% ==============================================================================
%> @brief Reflective boundary condition.
%
%> 
% ==============================================================================
classdef Reflective < BoundaryCondition
    % Base boundary condition class.
    %
    % This stores and computes the boundary flux along a given surface.
    
    properties (Access = protected)
        d_octants
    end
    
    methods
        
        % ======================================================================
        %> @brief Class constructor
        %>
        %> More detailed description of what the constructor does.
        %>
        %> @param boundary      Boundary flux class.
        %> @param side          Surface identifier.
        %>
        %> @return Instance of the Reflective class.
        % ======================================================================
        function obj = Reflective(boundary, input,  mesh, quadrature, side)
            obj = obj@BoundaryCondition(boundary, input,  mesh, quadrature, side);

            % If I'm a vertical side, a exhange mu's, e.g. I->II
            % If I'm a horizontal side, eta's echange (e.g. I->IV)
            
            % I have a side.  What are my incident octants? And what are
            % my exiting octants?
            if mesh.DIM == 1
                if side == Mesh.LEFT
                    obj.d_octants(1, 1) =  1; % I go in octant one, the right
                    obj.d_octants(1, 2) =  2; % and I reflect what goes left
                end
                if side == Mesh.RIGHT
                    obj.d_octants(1, 1) =  2; % I go in octant two, the left
                    obj.d_octants(1, 2) =  1; % and I reflect what goes right
                end
            elseif mesh.DIM == 2
                % In 2-D, we reflect over a whole half space.  That is, if we
                % are the left side, our incident condition points in all
                % directions such that (+mu, +/- eta), i.e. rightward and either
                % up or down.  That is, my (+,+)=I and (+,-)=IV comes from
                % (-,+)=II and (-,-)=III, respectively.
                if side == Mesh.LEFT
                    obj.d_octants(1, 1) =  1; % I go into the right, and
                    obj.d_octants(1, 2) =  2; % grab what goes left.
                    obj.d_octants(2, 1) =  4; % That is, I go into (I, IV) and
                    obj.d_octants(2, 2) =  3; % get what goes out of (II, III)
                elseif side == Mesh.RIGHT
                    obj.d_octants(1, 1) =  2; % I go into the left, and
                    obj.d_octants(1, 2) =  1; % grab what goes right.
                    obj.d_octants(2, 1) =  3; % That is, I go into (II, III) and                 
                    obj.d_octants(2, 2) =  4; % get what goes out of (I, IV)
                elseif side == Mesh.BOTTOM
                    obj.d_octants(1, 1) =  1; % I go to the top, which is
                    obj.d_octants(2, 1) =  2; % octants 1 and 2
                    obj.d_octants(1, 2) =  4; % and get reflection 
                    obj.d_octants(2, 2) =  3; % from 4 and 3, respectively  
                elseif side == Mesh.TOP
                    obj.d_octants(1, 1) =  3; % I go to the bottom, which is
                    obj.d_octants(2, 1) =  4; % octants 3 and 4
                    obj.d_octants(1, 2) =  2; % and get reflection 
                    obj.d_octants(2, 2) =  1; % from 2 and 1, respectively  
                else
                    error('Wrong side for 2-D Reflective')
                end
            else
                error('Only DIM = 1, 2  supported')
            end
            
            
            
        end
        
        % ======================================================================
        %> @brief Update the boundary flux.
        %
        %> Finish me.
        % ======================================================================
        function obj = update(obj)
            
            for o = 1:length(obj.d_octants(:, 1))
                
                o_in  = obj.d_octants(o, 1); % Put fluxes IN this octant
                o_out = obj.d_octants(o, 2); % Get fluxes OUT of this octant
                
                %for a = 1:number_angles_octant(obj.d_quadrature)
                      
                    if obj.d_side == Mesh.LEFT || obj.d_side == Mesh.RIGHT 
                        % Get the OUTGOING fluxes for this octant and angle
                        %f = get_psi_v(obj.d_boundary, o_out, a, Boundary.OUT);
                        f = get_psi_v_octant(obj.d_boundary, o_out, Boundary.OUT);
                        %f = f*0 + 1;
                        % Get the OUTGOING fluxes for this octant and angle
                        %set_psi_v(obj.d_boundary, o_in, a, f, Boundary.IN);
                        set_psi_v_octant(obj.d_boundary, o_in, f, Boundary.IN);
                        
                    elseif obj.d_side == Mesh.TOP || obj.d_side == Mesh.BOTTOM 
                        % Get the OUTGOING fluxes for this octant and angle
                        %f = get_psi_h(obj.d_boundary, o_out, a, Boundary.OUT);
                        f = get_psi_h_octant(obj.d_boundary, o_out, Boundary.OUT);
                        % Get the OUTGOING fluxes for this octant and angle
                        %set_psi_h(obj.d_boundary, o_in, a, f, Boundary.IN);
                        set_psi_h_octant(obj.d_boundary, o_in, f, Boundary.IN);
                    else
                        error(' 3-D not done yet')
                    end
                    
                %end
                
            end
            
            
        end
        
    end
    
end