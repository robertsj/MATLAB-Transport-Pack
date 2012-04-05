%> @file  SurfaceNode.m
%> @brief SurfaceNode class definition.
% ==============================================================================
%> @brief 
%>
%> Nearly verbatim translation to MATLAB from P. Romano's Python code.
% ==============================================================================
classdef SurfaceNode < handle

    properties
       surface
       positiveSense
    end
    
    methods (Access = public)
        
        % ======================================================================
        %> @brief Class constructor
        %> @return Instance of the SurfaceNode class.
        % ======================================================================
        function self = SurfaceNode(surface, positiveSense)
            
            self.surface = surface;
            self.positiveSense = positiveSense;
        end
        
        function b = contains(self, location)
            b = self.surface.positiveSense(location) == self.positiveSense; 
        end
        
    end % public methods
    
end
        