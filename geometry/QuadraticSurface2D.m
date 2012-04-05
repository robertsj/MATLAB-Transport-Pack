%> @file  QuadraticSurface2D.m
%> @brief QuadraticSurface2D class definition.
% ==============================================================================
%> @brief Defines an implicit surface of the form
%>        F(x,y) = Ax^2 + By^2 + Cxy + Dx + Ey + F
%
%> Nearly verbatim translation to MATLAB from P. Romano's Python code.
% ==============================================================================
classdef QuadraticSurface2D < Surface2D
    
    properties (Access = protected)
        A_
        B_
        C_
        D_
        E_
        F_
    end
    
    methods (Access = public)
        
        % ======================================================================
        %> @brief Class constructor
        %> @return Instance of the Surface2D class.
        % ======================================================================
        function self = QuadraticSurface2D()
            self.A_ = 0;
            self.B_ = 0;
            self.C_ = 0;
            self.D_ = 0;
            self.E_ = 0;
            self.F_ = 0;
        end
        
        % ======================================================================
        %> @brief Evaluate F(x, y) at a given point.
        % ======================================================================
        function f = sense(self, r)
            x = r.x;
            y = r.y;
            f = self.A_*x.^2 + self.B_*y.^2 + self.C_.*x.*y + ...
                self.D_.*x + self.E_.*y + self.F_;
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
            A = self.A_;
            B = self.B_;
            C = self.C_;
            D = self.D_;
            E = self.E_;
            F = self.F_;

            % Put F(x,y) = 0 in form of ax^2 + bx + c = 0
            q = y0 - m*x0;
            a = A + B*m^2 + C*m;
            b = 2*B*m*q + C*q + D + E*m;
            c = B*q^2 + E*q + F;
            
            % Determine intersection
            quad = b^2 - 4*a*c;
            if quad < 0
                
                r1 = 0;
                r2 = 0;
                
            elseif quad == 0
                
                x = -b/(2*a);
                y = m*(x - x0) + y0;
                r1.x = x; 
                r1.y = y;
                r2 = 0;
                
            elseif quad > 0
                
                x1 = (-b + sqrt(quad))/(2*a);
                y1 = m*(x1 - x0) + y0;
                x2 = (-b - sqrt(quad))/(2*a);
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
            
            