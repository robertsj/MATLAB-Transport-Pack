%> @file  BoundaryConditionMOC.m
%> @brief Boundary condition class definition.
% ==============================================================================
%> @brief Base MOC boundary condition class.
%
%> This class contains several methods specific to MOC boundaries.
% ==============================================================================
classdef BoundaryConditionMOC < BoundaryCondition

    methods
        
        % ======================================================================
        %> @brief Class constructor
        %>
        %> More detailed description of what the constructor does.
        %>
        %> @param boundary      Boundary flux class.
        %> @param side          Surface identifier.
        %>
        %> @return Instance of the Boundary class.
        % ======================================================================
        function this = BoundaryConditionMOC(boundary, mesh, quadrature, side)
            this = this@BoundaryCondition(boundary, mesh, quadrature, side); 
            
 
        end
% 
%         % ======================================================================
%         %> @brief Update the boundary flux.
%         % ======================================================================
%         this = update(this);
        
        
        
    end
    
    methods (Access = protected)
       
        
    end
 
end