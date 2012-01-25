%> @file  VacuumTrack.m
%> @brief Boundary condition class definition.
% ==============================================================================
%> @brief Vacuum MOC boundary condition class.
% ==============================================================================
classdef VacuumTrack < BoundaryConditionMOC
    
    methods
        
        % ======================================================================
        %> @brief Class constructor
        %>
        %> More detailed description of what the constructor does.
        %>
        %> @param boundary      Boundary flux class.
        %> @param side          Surface identifier.
        %>
        %> @return Instance of the VacuumTrack class.
        % ======================================================================
        function this = VacuumTrack(boundary, mesh, quadrature, side)
            this = this@BoundaryConditionMOC(boundary, mesh, quadrature, side); 
            
 
        end

        % ======================================================================
        %> @brief Update the boundary flux.
        % ======================================================================
        function this = update(this)
            
            % Do nothing.
            
        end
        
        
        
    end
    
    methods (Access = protected)
       
        
    end
 
end