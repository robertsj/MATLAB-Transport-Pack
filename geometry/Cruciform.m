%> @file  Cruciform.m
%> @brief Cruciform class definition.
% ==============================================================================
%> @brief Defines an implicit surface for a cruciform pin.
% ==============================================================================
classdef Cruciform < Surface2D
    
    properties (Access = protected)

    end
    
    methods (Access = public)
        
        % ======================================================================
        %> @brief Class constructor
        %> @return Instance of the Cruciform class.
        % ======================================================================
        function self = Cruciform()
            % nothing yet
        end
        
        % ======================================================================
        %> @brief Evaluate F(x, y) at a given point.
        %> @return The sense (i.e. in or out)
        % ======================================================================
        function f = sense(self, r)
            x = r.x;
            y = r.y;
            
            if (x <= 0.185850)
                
                yd = sqrt(.034540-x.^2) + 0.53340;
                
            elseif (x > 0.18570 && x <= 0.53340)
                
                yd = (0.53340 - sqrt(.1207560-(x-0.53340).^2));
                
            elseif  (x > 0.53340)
                
                yd = sqrt(.034540  - (x-0.53340).^2);
                
            end  

        end
        
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
        
        
        
    end % public methods
    
end
            
            