%> @file  ReflectiveTrack.m
%> @brief Boundary condition class definition.
% ==============================================================================
%> @brief Vacuum MOC boundary condition class.
% ==============================================================================
classdef ReflectiveTrack < BoundaryConditionMOC
    
    properties 
        d_index
        d_number_tracks
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
        %> @return Instance of the ReflectiveTrack class.
        % ======================================================================
        function this = ReflectiveTrack(boundary, mesh, quadrature, side)
            this = this@BoundaryConditionMOC(boundary, mesh, quadrature, side); 
            
            %this.d_index = boundary.d_side_index{side};
            this.d_number_tracks = length(boundary.d_side_index{side}(:,1));
            build_index(this);
        end

        % ======================================================================
        %> @brief Update the boundary flux.
        % ======================================================================
        function this = update(this)
            
            % Loop through all the tracks I own, grap the flux of each
            % track's feeder, and update my own track flux.
            for i = 1:this.d_number_tracks
                
                for p = 1:number_polar(this.d_quadrature)
                    
                % Get the feeder flux.
                psi = get_single_psi(this.d_boundary, this.d_index(i, 4), ...
                    this.d_index(i, 5), p, this.d_index(i, 6), 2);
                
                % Set my flux.
                set_single_psi(this.d_boundary, this.d_index(i, 1), ...
                    this.d_index(i, 2), p, this.d_index(i, 3), psi, 1);
                
                end
                
            end

            
        end
        
        
        
    end
    
    methods (Access = protected)
       
        % ======================================================================
        %> @brief Build simple table of what I grab and what I set.
        % ======================================================================
        function this = build_index(this)
            
            % Presize my index.
            this.d_index = zeros(this.d_number_tracks, 6);
            
            % Go through all tracks that leave me.
            for i = 1:this.d_number_tracks
                
                % These specify my tracks.
                o = this.d_boundary.d_side_index{this.d_side}(i, 1);
                a = this.d_boundary.d_side_index{this.d_side}(i, 2);
                t = this.d_boundary.d_side_index{this.d_side}(i, 3);
                
                % Get who feeds me.
                o_f = this.d_boundary.d_bindex{o, a}(t, 1);
                a_f = this.d_boundary.d_bindex{o, a}(t, 2);
                t_f = this.d_boundary.d_bindex{o, a}(t, 3);
                
                % Add to the table.
                this.d_index(i, :) = [o a t o_f a_f t_f];
                
            end
        end
        
    end
 
end