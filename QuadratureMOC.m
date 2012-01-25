%> @file  QuadratureMOC.m
%> @brief QuadratureMOC class definition.
% ===========================================================================
%> @brief Base quadrature class for MOC.
%
%> For the polar angle, T-Y quadrature for orders 1, 2, and 3 is allowed.
%>
%> Azimuthal angles are defined by derived classes.
%> 
%> The azimuthal angles are stored for the first quadrant.  Other
%> quadrants are given an additive constant.
%>
% ===========================================================================
classdef QuadratureMOC < handle
    
    properties
        %> Azimuthal angles (positive quadrant only)
        d_phi
        %> Azimuthal weights
        d_weight_phi
        %> Polar cosines (w/r to the x-y plane)
        d_mu
        %> Polar weights
        d_weight_mu
        %> Number of azimuthal angles
        d_number_azimuth
        %> Number of polar angles
        d_number_polar
        %> Number horizontal intercepts
        d_number_x
        %> Number vertical intercepts
        d_number_y
        %>
        d_number_tracks
        %> Number of tracks
        d_number_space
        %>
        d_enter
        d_exit
        d_space
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
    
    methods
        
        % ======================================================================
        %> @brief Class constructor
        %
        %> @param number_polar  Number of polar angles used.
        %> @return              Instance of the QuadratureMOC class.
        % ======================================================================
        function obj = QuadratureMOC(number_polar)
            % Tabuchi-Yamamoto polar quadrature.
            if number_polar == 1
                obj.d_mu = 0.798184;
                obj.d_weight_mu = 2.0;
                
            elseif number_polar == 2
                obj.d_mu = [0.363900; 0.899900];
                obj.d_weight_mu = 2*[0.212854; 0.787146];
                
            elseif number_polar == 3
                obj.d_mu = [0.166648; 0.537707; 0.932954];
                obj.d_weight_mu = 2*[0.046233; 0.283619; 0.670148];
            else
                error('Invalid polar order');
            end
            obj.d_number_polar = number_polar;
        end
        
        function p = phi(obj, o, i)
            p = obj.d_phi(i) + (o-1)*pi/2;
        end
        
        function w = weight_phi(obj, i)
            w = obj.d_weight_phi(i); 
        end
        
        function m = mu(obj, i)
            m = obj.d_mu(i); 
        end
        
        function w = weight_mu(obj, i)
            w = obj.d_weight_mu(i); 
        end
        
        function n = number_azimuth(obj)
            n = obj.d_number_azimuth * 4;
        end
        
        function n = number_angles_octant(obj)
        	n = obj.d_number_azimuth;
        end
        
        function n = number_angles(obj)
            n = obj.d_number_azimuth;
        end
        
        function n = number_polar(obj)
            n = obj.d_number_polar;
        end
        
        function n = number_x(obj, m)
            n = obj.d_number_x(m);
        end
        
        function n = number_y(obj, m)
            n = obj.d_number_y(m);
        end
        
        function n = number_tracks(obj, m)
            n = obj.d_number_tracks(m);
        end     
        
        function n = total_number_track(obj)
            n = sum(obj.d_number_tracks(:));
        end
        
        % ======================================================================
        %> @brief Computes cardinal angle index.
        %> @return  Index.
        % ======================================================================
        function y = index(obj, o, a)
            y = a + (o - 1) * number_angles_octant(obj);
        end
        
        function p = uniform(obj, a, b, m, f)
            pp = linspace(a, b, m+1);
            p = pp(1:end-1)+pp(2)/2;
            if f == 1
                p = p(end:-1:1);
            end
        end
        
    end
    
end