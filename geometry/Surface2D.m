%> @file  Surface2D.m
%> @brief Surface2D class definition.
% ==============================================================================
%> @brief Defines an implicit surface.
%
%> Nearly verbatim translation to MATLAB from P. Romano's Python code.
% ==============================================================================
classdef Surface2D < handle
    
    properties (Access = protected)

    end
    
    methods (Access = public)
        
%         % ======================================================================
%         %> @brief Class constructor
%         %> @return Instance of the Surface2D class.
%         % ======================================================================
%         function self = Surface2D()
%             % nothing yet
%         end
        
        % ======================================================================
        %> @brief Evaluate F(x, y) at a given point.
        %> @return The sense (i.e. in or out)
        % ======================================================================
        f = sense(self, r)
        
        % ======================================================================
        %> @brief Solves the equation F(r0 + s*omega) = 0 for the point s
        %> 
        %> This is the point of intersection of the surface F(x,y) and
        %> the line r = r0 + s*Omega.
        % ======================================================================
%         function d = distance(self, r, omega)
%             
%         end
        
        % ======================================================================
        %> @brief Determine if location r is 'above' or 'below' the surface
        %>        by determining the sense
        % ======================================================================
        function b = positiveSense(self, r)  
            b = self.sense(r) > 0; 
        end
        
        % ======================================================================
        %> @brief Determine the points of intersection of the surface with a
        %>        line drawn defined from point r with angle phi. 
        %>
        %>  Effectively solves the equation F(x, y0+m*(x-x0)) = 0 for x.
        % ======================================================================
        [r1, r2] = intersectLine(self, r, phi)
        
        
        % ======================================================================
        %> @brief Draw me.
        % ======================================================================
        function draw(self, bound, n)
            if nargin == 2
                n = 100;
            end
            r.x = rand(n, 1)*(bound(2)-bound(1))+bound(1);
            r.y = rand(n, 1)*(bound(2)-bound(1))+bound(1);
            b = self.positiveSense(r);
            x_i = r.x(find(~b));
            y_i = r.y(find(~b));
            x_o = r.x(find(b));
            y_o = r.y(find(b));
            plot(x_i,y_i,'b.', x_o, y_o, 'k.')
        end
        
    end % public methods
    
end
            
            