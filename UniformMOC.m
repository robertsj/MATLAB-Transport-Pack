%> @file  UniformMOC.m
%> @brief UniformMOC class definition.
% ===========================================================================
%> @brief Uniform azimuthal quadrature.
%
%> For the polar angle, T-Y quadrature for orders 1, 2, and 3 is allowed.
%>
%> Azimuthal angles are defined by derived classes.
%> 
%> The azimuthal angles are stored for the first quadrant.  Other
%> quadrants are given an additive constant.
%>
% ===========================================================================
classdef UniformMOC < QuadratureMOC
    
    properties

    end
    
    
    methods
        
        function obj = UniformMOC(number_azimuth, number_polar, number_space)
            
            % Call base class
            obj = obj@QuadratureMOC(number_polar);

            obj.d_number_azimuth = number_azimuth;
            obj.d_number_space   = number_space;
            obj.d_phi      = zeros(number_azimuth, 1);
            obj.d_number_x = zeros(number_azimuth, 1);
            obj.d_number_y = zeros(number_azimuth, 1);
            
            % Four way angular symmetry
            delta = 0.5 * pi / number_azimuth;
            
            for m = 1:number_azimuth
                
                % Preliminary azimuth
                if ( m <= number_azimuth * 0.5 )
                    phi = delta * (m - 0.5);
                else
                    %error(' oops ')
                    phi = delta * (number_azimuth + 0.5 -m);
                end
                % Number of x and y intercepts
                tan_phi  = tan(phi);
                number_x = abs(ceil(number_space * tan_phi / (tan_phi + 1.0)));
                number_y = number_space - number_x;
                
                % Actual Azimuth
                if (m <= number_azimuth * 0.5)
                    phi = atan(number_x / number_y);
                elseif (m > number_azimuth * 0.5) % 2nd quadrant
                    phi = pi/2 - atan(number_x / number_y);
                else
                    error (' Error in computing actual azimuths... ');
                end
                obj.d_phi(m) = phi;
            end
            % Calculating azimuthal weights
            obj.d_weight_phi = zeros(number_azimuth, 1);
            obj.d_weight_phi(1) = (0.5*pi + obj.d_phi(2) - ...
                obj.d_phi(number_azimuth) ) / pi;
            obj.d_weight_phi(number_azimuth) = ...
                (0.5*pi + obj.d_phi(1) - obj.d_phi(number_azimuth-1)) / pi;
            % All other m
            for m = 2:(number_azimuth-1)
                obj.d_weight_phi(m) = (obj.d_phi(m+1) - obj.d_phi(m-1)) / pi;
            end
        end
        
        function p = phi(obj, o, i)
            p = obj.d_phi(i) + (o-1)*pi/2;
        end
        
        function m = mu(obj, i)
            m = obj.d_mu(i); 
        end
        
    end
    
end