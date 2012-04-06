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
        function this = VacuumTrack(boundary, input, mesh, quadrature, side)
            this = this@BoundaryConditionMOC(boundary, input, mesh, quadrature, side); 
            
 
        end
        
        % ======================================================================
        %> @brief Set the boundary flux.
        %
        %> Note, nothing needs to be done for reflective, since the initial
        %> guess is zero.
        % ======================================================================
        function this = set(this)
            % do nothing
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