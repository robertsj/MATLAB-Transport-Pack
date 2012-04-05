%> @file  Geometry.m
%> @brief Geometry class definition.
% ==============================================================================
%> @brief The Geometry class contains information on the domain, geometry,
%>        regions, and materials in the problem.
%
%> Nearly verbatim translation to MATLAB from P. Romano's Python code.
% ==============================================================================
classdef Geometry < handle
    
    properties (Access = protected)
        %> Width of the problem domain
        width
        %> Height of the problem domain
        height
        %> Cell array of @ref Region objects in the geometry.
        regions
        %> Map of regions by name
        regionsByName
    end
    
    methods (Access = public)
        
        % ======================================================================
        %> @brief Class constructor
        %> @return Instance of the Geometry class.
        % ======================================================================
        function self = Geometry(width, height)
            self.width         = width;
            self.height        = height;
            self.regions       = cell(0);
            self.regionsByName = containers.Map();
        end
        
        
        % ======================================================================
        %> @brief Add a @ref Region object to the geometry, storing it in the
        %>        regions list and regionsByName dictionary.
        % ======================================================================
        function self = addRegion(self, region)
            self.regions{length(self.regions)+1} = region;
            self.regionsByName(region.name) = region;
        end
            
        % ======================================================================
        %> @brief For a given x,y coordinate and angle, determine the
        %>        coordinates of the outgoing boundary point
        % ======================================================================
        function r = endpoint(self, r, phi)
            
            w = self.width;
            h = self.height;
            x = r.x;
            y = r.y;
            m = tan(phi);
            
            % Determine points of intersection with boundaries.
            points = [ 0,            y - x*m
                       w,            y + (w-x)*m
                       x - y/m,      0
                       x + (h-y)/m,  h           ];  
            
            % Determine which point is on boundary.
            for i = 1:length(points(:, 1))
                % Skip the trivial pair (x,y).
                if points(i, 1) == x && points(i, 2) == y
                    continue
                end
                if points(i, 1) >= 0 && points(i, 1) <= w && ...
                   points(i, 2) >= 0 && points(i, 2) <= h 
                   r.x = points(i, 1);
                   r.y = points(i, 2);
                end
            end
            
        end
        
        %> @brief Sketch the geometry.
        function sketch(self, n)
            x = rand(n, 1)*self.width;
            y = rand(n, 1)*self.height;
            ridx = zeros(n, 1);
%             for i = 1:n
%                 r.x = x(i);
%                 r.y = y(i);
%                 found = 0;
%                 for reg = 1:length(self.regions)
%                     if self.regions{reg}.node.contains(r)
%                         ridx(i) = reg;
%                         found = 1;
%                         break;
%                     end
%                 end
%                 if ~found
%                     disp(['uncontained point, x = ',num2str(r.x), ' y = ', num2str(r.y)])
%                 end
%             end
            r.x = x;
            r.y = y;
            colors = rand(length(self.regions), 3);
            clf
            hold on
            for reg = 1:length(self.regions)
                b = self.regions{reg}.node.contains(r);
                xx = x(find(b));
                yy = y(find(b));
                plot(xx, yy, '.', 'Color', [colors(reg,:)], 'MarkerSize',12);
            end
            
            
%             clf
%             hold on
%             for i = 1:n
%                 if ridx(i) > 0
%                 plot(x(i), y(i), '.', 'Color', [colors(ridx(i),:)], 'MarkerSize',12);
%                 end
%             end
            axis([0 self.width 0 self.height])
            axis square
        end
        
    end % public methods
    
end
            
            