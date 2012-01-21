%> @file  UniformMOC.m
%> @brief UniformMOC class definition.
% ===========================================================================
%> @brief Uniform azimuthal quadrature.
%
%> This quadrature takes a user-defined number of tracks, used for all
%> angles, and a number of azimuthal angles per quadrant.  The azimuthal
%> angles are initial uniformly distributed in the first quadrant.  This
%> is used to select the number of x and y intercepts for the given
%> angle.  These intercepts are distributed evenly on the surfaces.  The
%> azimuthal angles are then adjusted to match the x and y intercepts.
%> Note, this means some combinations do not work.  If an error about
%> repeated angles arises, the user must choose fewer angles or more
%> spatial points.  Generally, for an even number of angles and tracks,
%> the spacing
%> should be about 3 times higher to ensure they don't repeat.  For an odd
%> number of angles and track, a factor of 2 plus 1 seems to work.
%>
% ===========================================================================
classdef UniformMOC < QuadratureMOC
    
    properties

    end
    
    
    methods
        
        function obj = UniformMOC(number_azimuth, number_polar, number_space)
            
            % For now, require odd number of phi's
            
            
            % Call base class
            obj = obj@QuadratureMOC(number_polar);

            obj.d_number_azimuth = number_azimuth;
            obj.d_number_space   = number_space;
            
            % Define evenly spaced azimuthal angles over [0, pi/2].
            a = pi/4/number_azimuth;
            b = pi/2-pi/4/number_azimuth;
            obj.d_phi = linspace(a, b, number_azimuth);
            
            obj.d_number_x      = zeros(number_azimuth, 1);
            obj.d_number_y      = zeros(number_azimuth, 1);
            obj.d_number_tracks = zeros(number_azimuth, 1);
            
            % Four way angular symmetry
            delta = 0.5 * pi / number_azimuth;
            
            for m = 1:number_azimuth

                % Number of x and y intercepts
                tan_phi  = tan(obj.d_phi(m));
                % We do this if statement so that the adjusted angles
                % remain symmetric about pi/4.  Asthetics.
                if m <= number_azimuth/2
                    number_x = abs(ceil(number_space * tan_phi ...
                        / (tan_phi + 1.0)));
                else
                    number_x = abs(floor(number_space * tan_phi ...
                        / (tan_phi + 1.0)));  
                end
                number_y = number_space - number_x;
                
                % Actual Azimuth
                phi = atan(number_x / number_y);
                obj.d_phi(m) = phi;
                obj.d_number_x(m) = number_x;
                obj.d_number_y(m) = number_y;
                obj.d_number_tracks(m) = number_x + number_y;
            end
            
            % For now, don't do any a priori checks.  Basically, need
            % enough spatial points to allow the angles wiggle room.
            if length(unique(obj.d_phi)) ~= length(obj.d_phi)
                error('Repeat angles!')
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
            
            % Calculate intercepts on a square.
            obj.d_enter = cell(number_azimuth, 1);
            obj.d_exit  = cell(number_azimuth, 1);
            
            for m = 1:number_azimuth
                % Number of tracks
                nx = obj.d_number_x(m);
                ny = obj.d_number_y(m);
                n  = nx+ny;
                % Horizontal and vertical steps.
                dx = 1 / nx;
                dy = 1 / ny;
                % Perpendicular distance between tracks
                obj.d_space(m) = 1 / cos(obj.d_phi(m));
                % First quadrant only (0, pi/2)
                %   Uniformly spaced entrances
                enters(1:ny, 1)   = 0.0;
                enters(1:ny, 2)   = (1-dy/2):(-dy):(0);
                enters(ny+1:n, 1) = (dx/2):(dx):(1);
                enters(ny+1:n, 2) = 0.0;   
                %   The exits are similar.
                exits(1:nx, 1) = (dx/2):(dx):(1);
                exits(1:nx, 2) = 1;
                exits(nx+1:n, 1) = 1;
                exits(nx+1:n, 2) = (1-dy/2):(-dy):(0); 
                obj.d_enter{m} = enters;
                obj.d_exit{m}  = exits;
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