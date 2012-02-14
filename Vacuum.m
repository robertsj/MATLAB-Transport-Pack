%> @file  Vacuum.m
%> @brief Vacuum boundary condition class definition.
% ==============================================================================
%> @brief Vacuum boundary condition.
%
%> This is essentially an empty condition, since nothing needs to be done.
% ==============================================================================
classdef Vacuum < BoundaryCondition
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
        %> @return Instance of the Vacuum class.
        % ======================================================================
        function obj = Vacuum(boundary, input, mesh, quadrature, side)
            obj = obj@BoundaryCondition(boundary, input, mesh, quadrature, side); 

        end
        
        % ======================================================================
        %> @brief Update the boundary flux.
        %
        %> Note, nothing needs to be done for vacuum, as the incident fluxes are
        %> initialized as zero.
        % ======================================================================
        function obj = update(obj)
            % do nothing
        end
        
    end
    
end