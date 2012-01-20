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
        %> Number of tracks
        d_number_space
    end
    
    
    methods
        
        function obj = QuadratureMOC(number_polar)
            
            if number_polar == 1
                obj.d_mu = 0.798184;
                obj.d_weight_mu = 1.0;
                
            elseif number_polar == 2
                obj.d_mu = [0.363900; 0.899900];
                obj.d_weight_mu = [0.212854; 0.787146];
                
            elseif number_polar == 3
                obj.d_mu = [0.166648; 0.537707; 0.932954];
                obj.d_weight_mu = [0.046233; 0.283619; 0.670148];
            else
                error('Invalid polar order');
            end
            obj.d_number_polar = number_polar;
        end
        
        function p = phi(obj, o, i)
            p = obj.d_phi(i) + (o-1)*pi/2;
        end
        
        function m = mu(obj, i)
            m = obj.d_mu(i); 
        end
        
    end
    
end