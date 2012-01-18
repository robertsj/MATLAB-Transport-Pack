%> @file  Reflective.m
%> @brief Reflective boundary condition class definition.
% ==============================================================================
%> @brief Reflective boundary condition.
%
%> This is essentially an empty condition, since nothing needs to be done.
% ==============================================================================
classdef Reflective < BoundaryCondition
    % Base boundary condition class.
    %
    % This stores and computes the boundary flux along a given surface.
    
    properties (Access = protected)
        
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
        function obj = Reflective(boundary, side)
            obj = obj@Boundary(boundary, side); 

        end
        
        % ======================================================================
        %> @brief Update the boundary flux.
        %
        %> Finish me.
        % ======================================================================
        function obj = update(obj)
            % do nothing
        end
        
    end
    
end