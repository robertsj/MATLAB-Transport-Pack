%> @file  Quadrature.m
%> @brief Quadrature class definition.
% ===========================================================================
%> @brief Base quadrature class.
%
%> And here we can put some more detailed informations about the class.
%> All quadratures must be ordered so that the signs of the cosines are 
%> arranged in the following manner:
%>
%>    indices | mu  | eta | xi  
%>    ----------------------------------
%>       1: N | +   |  +  | +        (first octant)
%>     N+1:2N | -   |  +  | +        (second octant)
%>    2N+1:3N | -   |  -  | +        (third octant)
%>    3N+1:4N | +   |  -  | +        (fourth octant)
%>    4N+1:5N | +   |  +  | -        ...
%>    5N+1:6N | -   |  +  | -       
%>    6N+1:7N | -   |  -  | -      
%>    7N+1:8N | +   |  -  | -       
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
        octant = [  1  1  1 
                   -1  1  1  
                   -1 -1  1
                    1 -1  1
                    1  1 -1
                   -1  1 -1
                   -1 -1 -1
                    1 -1 -1  ];
    end
    
    properties (Access = protected)
        %> Quadrature weights
        d_weights 
        %> x-axis cosines
        d_mu
        %> y-axis cosines
        d_eta
        %> z-axis cosines
        d_xi
        %> type of quadrature being used
        d_name
        %> number of angles
        d_number_angles
        %> number of octants
        d_number_octants
        %> problem dimension
        d_dim
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
        function obj = Quadrature(order, dim)
            if (order < 1)
                error('Quadrature order must be greater than 1.')
            end
            obj.d_order = order;
            obj.d_dim = dim;
            if dim == 1
                obj.d_number_octants = 2;
            elseif dim == 2
                obj.d_number_octants = 4;
            elseif dim == 3
                obj.d_number_octants = 8;
            else
                error('Invalid dimension.')
            end
                
        end
        
    end
    
    methods (Access = public)
        
        % ----------------------------------------------------------------------
        % Public Interface
        % ----------------------------------------------------------------------
        
        % ======================================================================
        %> @brief Return weights for all angles.
        %> @return  Vector of weights.
        % ======================================================================
        function w = weights(obj)
            w = obj.d_weights;
        end
        
        % ======================================================================
        %> @brief Return weight for an angle (indexed by cardinal index)
        %> @return  Vector of weights.
        % ======================================================================
        function w = weight(obj, a)
            w = obj.d_weights(a);
        end
        
        % ======================================================================
        %> @brief Return all angles.
        %> @return  Vector of weights.
        % ======================================================================
        function [mu, eta, xi] = angles(obj)
            mu  = obj.d_mu;
            eta = obj.d_eta;
            xi  = obj.d_xi;
        end
        
        % ======================================================================
        %> @brief Return angles for a given cardinal index.
        %> @return  Set of angles.
        % ======================================================================
        function [mu, eta, xi] = angle(obj, o, a)
            check_index(obj, a);
            mu  = obj.d_mu(a)  * obj.octant(o, 1);
            if obj.d_dim > 1
            	eta = obj.d_eta(a) * obj.octant(o, 2);
            else
                eta = 0;
            end
            if obj.d_dim > 2
            	xi = obj.d_eta(a) * obj.octant(o, 3);
            else
                xi = 0;
            end
        end   
        
        % ======================================================================
        %> @brief Total number of angles.
        %> @return  Number of angles.
        % ======================================================================
        function y = number_angles(obj)
            y = obj.d_number_angles;
        end
        
        % ======================================================================
        %> @brief Total number of angles in an octant.
        %> @return  Number of angles.
        % ======================================================================
        function y = number_angles_octant(obj)
            y = obj.d_number_angles/obj.d_number_octants;
        end     
        
        % ======================================================================
        %> @brief Total number of octants.
        %> @return  Number of octants.
        % ======================================================================
        function y = number_octants(obj)
            y = obj.d_number_octants;
        end   
        
        % ======================================================================
        %> @brief Angle bounds within an octant.
        %> @return  Bounds.
        % ======================================================================
        function b = bounds(obj, o)
            DBC.Require('o > 0 && o <=  number_octants(obj)') 
            n_a = number_angles_octant(obj);
            b(1) = (o - 1)*n_a + 1; 
            b(2) = n_a * o; 
        end
        
        % ======================================================================
        %> @brief Computes cardinal angle index.
        %> @return  Index.
        % ======================================================================
        function y = index(obj, o, a)
            y = a + (o - 1) * number_angles_octant(obj);
        end
        
    end
    
    methods (Static)
        
        function y = angular_norm(dim)
            y = 1.0 / 4.0 / pi;
            if dim == 1
                y = 1.0 / 2.0;
            end
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



