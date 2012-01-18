%> @file  Quadrature.m
%> @brief Quadrature class definition.
% ===========================================================================
%> @brief Base quadrature class.
%
%> And here we can put some more detailed informations about the class.
%> All quadratures must be ordered so that the signs of the cosines are 
%> arranged in the following manner:
%>
%>    indices | mu(range) | eta(indices)
%>    ----------------------------------
%>       1: N |    +      |      +           (first octant)
%>     N+1:2N |    -      |      +           (second octant)
%>    2N+1:3N |    -      |      -           (third octant)
%>    3N+1:4N |    +      |      -           (fourth octant)
%>
%> Note that N is the number of angles per quadrant.  Outside of the given
%> pattern, the angles need only be consistently ordered, i.e.
%> 
%>   abs(mu(i*N+1)) = abs(mu(j*N+1)) for i,j = 0, 1, 2, 3
%>
%> though decreasing absolute value is suggested.
% ===========================================================================
classdef Quadrature
    
    properties
        d_order
    end
  
    properties (Constant)
        %> Octant cosign signs.  
        octant = [1 1; -1 1; -1 -1; 1 -1];      
        %> Number of octants.
        number_octants = 4;
    end
    
    properties (Access = protected)
        %> Quadrature weights
        d_weights 
        %> x-axis cosines
        d_mu
        %> y-axis cosines
        d_eta
        %> type of quadrature being used
        d_name
        %> number of angles
        d_number_angles
    end
    
    methods
        
        % ======================================================================
        %> @brief Class constructor
        %>
        %> More detailed description of what the constructor does.
        %>
        %> @param order         Quadrature order. This differs from quadrature 
        %>                      to quadrature, but e.g. for level symmetric, 
        %>                      it's the number of unique directional cosines.
        %>
        %> @return Instance of the Quadrature class.
        % ======================================================================
        function obj = Quadrature(order)
            if (order < 1)
                error('Quadrature order must be greater than 1.')
            end
            obj.d_order = order;
        end
        
    end
    
    methods (Access = public)
        
        % ----------------------------------------------------------------------
        % Public Interface
        % ----------------------------------------------------------------------
        
        function w = weights(obj)
        % function w = weights(obj)
        %   Returns the quadrature weights.
            w = obj.d_weights;
        end
        
        function w = weight(obj, a)
        % function w = weight(obj)
        %   Returns the quadrature weight for a given cardinal index a.  
            %check_index(obj, a);
            w = obj.d_weights(a);
        end
        
        function [mu, eta] = angles(obj)
        % function [mu, eta] = angles(obj);
        %   Returns the quadrature points.
            mu = obj.d_mu;
            eta = obj.d_eta;
        end
        
        function [mu, eta] = angle(obj, o, a)
        % function [mu, eta] = angle(obj, a);
        %   Returns the quadrature points for a given cardinal index a
            check_index(obj, a);        
            mu  = obj.d_mu(a)  * obj.octant(o, 1);
            eta = obj.d_eta(a) * obj.octant(o, 2);
        end   
        
        function y = number_angles(obj)
        % function y = number_angles(obj)
        %   Return the number of angles in the quadrature.
            y = obj.d_number_angles;
        end
        
        function y = number_angles_octant(obj)
        % function y = number_angles(obj)
        %   Return the number of angles in the quadrature.
            y = obj.d_number_angles/4;
        end      
        
        
        function b = bounds(obj, o)
            n = number_angles_octant(obj);
            if o == 1
                b(1) = 1;
                b(2) = n;
            elseif o == 2
                b(1) = n + 1;
                b(2) = 2*n;
            elseif o == 3
                b(1) = 2*n + 1;
                b(2) = 3*n;
            else
                b(1) = 3*n + 1;
                b(2) = 4*n;
            end
        end
        
        function y = index(obj, o, a)
            y = a + (o - 1) * number_angles_octant(obj);
        end
        
    end
    
    methods(Access = protected)
        function bool = check_index(obj, a)
        % function bool = check_index(obj, angle);
        %   Checks that angle is a value index into the quadrature
            if a > length(obj.weights) || a < 1
                error('Invalid angle index.')
            else
                bool = 1;
            end
        end 
    end
end



