%> @file  Plane.m
%> @brief Plane class definition.
% ==============================================================================
%> @brief Defines an implicit surface of the form
%>        F(x,y) = Ax + By + C = 0
%
%> Nearly verbatim translation to MATLAB from P. Romano's Python code.
% ==============================================================================
classdef Plane < QuadraticSurface2D

    methods (Access = public)
        
        % ======================================================================
        %> @brief Class constructor
        %> @return Instance of the Surface2D class.
        % ======================================================================
        function self = Plane(xc, yc, const)
            self = self@QuadraticSurface2D();
            self.D_ = xc;
            self.E_ = yc;
            self.F_ = const;
        end
   
    end % public methods
    
end
        