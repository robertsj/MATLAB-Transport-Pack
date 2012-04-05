%> @file  Circle.m
%> @brief Plane class definition.
% ==============================================================================
%> @brief Defines an implicit surface of the form
%>        F(x,y) = (x-x0)^2 + (y-y0)^2 - R^2 = 0
%
%> Nearly verbatim translation to MATLAB from P. Romano's Python code.
% ==============================================================================
classdef Circle < QuadraticSurface2D
    
    methods (Access = public)
        
        % ======================================================================
        %> @brief Class constructor
        %> @return Instance of the Circle class.
        % ======================================================================
        function self = Circle(x, y, R)
            self.A_ = 1;
            self.B_ = 1;
            self.C_ = 0;
            self.D_ = -2*x;
            self.E_ = -2*y;
            self.F_ = x^2 + y^2 - R^2;
        end
        
        function v = x(self)
            v = -0.5 * self.D_;
        end
        
        function v = y(self)
            v = -0.5 * self.E_;
        end
        
        function v = R(self)
            v = sqrt(self.x()^2 + self.y()^2 - self.F_);
        end
        
        function a = area(self)
            a = pi * self.R()^2; 
        end
        
        % ======================================================================
        %> @brief Determine the points of intersection of the surface with a
        %>        line drawn defined from point r with angle phi. 
        %>
        %>  Effectively solves the equation F(x, y0+m*(x-x0)) = 0 for x.
        % ======================================================================
        function [r1, r2] = intersectLine(self, r, phi)
            
            % Reference point and angle for line
            x0 = r.x;
            y0 = r.y;
            m = tan(phi);

            % Center of circle
            xc = self.x();
            yc = self.y();
            R  = self.R();

            % Put F(x,y) = 0 in form of ax^2 + bx + c = 0
            a = 1 + m^2;
            k = -xc - m^2*x0 + m*(y0 - yc);
            c = xc^2 + (m*x0)^2 - 2*m*x0*(y0 - yc) + (y0 - yc)^2 - R^2;
            
            % Determine intersection
            quad = k^2 - a*c;
            if quad < 0
                
                r1 = 0;
                r2 = 0;
                
            elseif quad == 0
                
                x = -k/a;
                y = m*(x - x0) + y0;
                r1.x = x; 
                r1.y = y;
                r2 = 0;
                
            elseif quad > 0
                
                x1 = (-k + sqrt(quad))/a;
                y1 = m*(x1 - x0) + y0;
                x2 = (-k - sqrt(quad))/a;
                y2 = m*(x2 - x0) + y0;
                d1 = abs(y1-y0);
                d2 = abs(y2-y0);
                if d1 < d2
                    r1.x = x1; r1.y = y1;
                    r2.x = x2; r2.y = y2;
                else
                    r1.x = x2; r1.y = y2;
                    r2.x = x1; r2.y = y1;
                end
                
            end
            
        end
        
    end % public methods
    
end
