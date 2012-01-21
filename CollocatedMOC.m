%> @file  CollocatedMOC.m
%> @brief CollocatedMOC class definition.
% ===========================================================================
%> @brief Collocated azimuthal quadrature.
%
%> This quadrature has a number of unique features.  It satisfies cyclic
%> tracking requirements in a rectangular region.  It also uses a finite 
%> set of points along a surface as track entrance (and exit) points. In
%> first implementation, we consider only square regions.
%>
%> The user sets the number of spatial points along a side.  Currently,
%> this must be a power of 3, and only 9 and 27 are supported.  
%> 
%> The angles (between \f$ \pi/4 \f$ and \f$ \pi/2 \f$) are uniquely defined 
%> by
%> \f[
%>      \tan{\phi_i} = 3^i \, , \,\,\, i = 0, \, \cdots n \, ,
%> \f] 
%> where 
%> \f[
%>      n = \mathrm{floor}(\log_3(N)+1) \,
%> \f]
%> and  \em N is the number of spatial points. Note that n/N is
%> maximized when \em N  is a power of three.  Using 9 spatial points
%> surface yields 5 angles per quadrant.  27 points yields 7 angles per
%> quadrant.  The other angles are defined by symmetry about \f$ \pi/4 \f$.
%>
%> For weights, the user can use an arc length approximation or a 
%> quadrature rule defined for the points that may provide better 
%> accuracy.
% ===========================================================================
classdef CollocatedMOC < QuadratureMOC
    
    properties

    end
    
    
    methods
        
        function obj = CollocatedMOC(number_space, ...
                                     number_polar, ...
                                     weight_flag)
            
            % For now, require odd number of phi's
            
            
            % Call base class
            obj = obj@QuadratureMOC(number_polar);
            obj.d_number_space   = number_space;
            
            if (number_space == 3)
                obj.d_number_azimuth = 3;            
            elseif (number_space == 9)
                obj.d_number_azimuth = 5;
            elseif (number_space == 27)
                obj.d_number_azimuth = 7;
            elseif (number_space == 81)
                obj.d_number_azimuth = 9;                
            else
                error('invalid space number')
            end
            
            obj.d_number_x      = number_space*ones(obj.d_number_azimuth, 1);
            obj.d_number_y      = number_space*ones(obj.d_number_azimuth, 1);
            obj.d_number_tracks = zeros(obj.d_number_azimuth, 1);
            
            obj.d_phi = zeros(obj.d_number_azimuth, 1);
            obj.d_phi(ceil(obj.d_number_azimuth/2)) = pi/4;
            for i = 1:floor(obj.d_number_azimuth/2)
                obj.d_phi(ceil(obj.d_number_azimuth/2)+i) = atan(3^i);
                obj.d_number_y(obj.d_number_azimuth-i+1) = 3^(i-1);
                obj.d_number_x(i) = 3^(i-1);
                obj.d_phi(ceil(obj.d_number_azimuth/2-i)) = atan(1/3^i);
                
            end
            
            % Compute weights.
            
            if weight_flag == 0
                
                u = zeros(obj.d_number_azimuth, 1);
                A = zeros(obj.d_number_azimuth);
                p = obj.d_phi*2/pi; % scale phi to be in (0, 1)
                for i = 1:obj.d_number_azimuth
                    u(i) = 1/i; % integrated monomial over (0, 1)
                    for j = 1:obj.d_number_azimuth
                        A(i, j) = p(j)^(i-1);
                    end
                end
                obj.d_weight_phi = A\u;
                
            else
                % Calculating azimuthal weights
                obj.d_weight_phi    = zeros(obj.d_number_azimuth, 1);
                obj.d_weight_phi(1) = (0.5*pi + obj.d_phi(2) - ...
                    obj.d_phi(obj.d_number_azimuth) ) / pi;
                obj.d_weight_phi(obj.d_number_azimuth) = ...
                    (0.5*pi + obj.d_phi(1) - ...
                    obj.d_phi(obj.d_number_azimuth-1)) / pi;
                % All other m
                for m = 2:(obj.d_number_azimuth-1)
                    obj.d_weight_phi(m) = (obj.d_phi(m+1) - ...
                        obj.d_phi(m-1)) / pi;
                end
            end

            % Reference values
            %
            %   quad.7 = [0.065371507540233 0.030298024005786   
            %             0.243764494414337 0.321131948079288
            %   quad.5 = [0.145786094758057 0.169492101633954 
            %             0.369443607215977]
            %
            %   arc.7  = [0.047007156347619 0.090632513479225 ...
            %             0.214776712522723 0.295167235300867]
            %   arc.5  = [0.137639669826844 0.214776712522723 
            %             0.295167235300867]
            %
            
            % Calculate intercepts on a square.
            obj.d_enter = cell(obj.d_number_azimuth, 1);
            obj.d_exit  = cell(obj.d_number_azimuth, 1);            
            for m = 1:obj.d_number_azimuth
                % Number of tracks
                nx = obj.d_number_x(m);
                ny = obj.d_number_y(m);
                n  = nx+ny;
                obj.d_number_tracks(m) = n;
                % Horizontal and vertical steps.
                dx = 1 / nx;
                dy = 1 / ny;
                % Perpendicular distance between tracks
                obj.d_space(m) = 1 / cos(obj.d_phi(m));
                % First quadrant only (0, pi/2)
                %   Uniformly spaced entrances
                enters(1:ny, 1)   = 0.0;
                enters(1:ny, 2)   = uniform(obj, 0, 1, ny, 1);
                enters(ny+1:n, 1) = uniform(obj, 0, 1, nx, 0);
                enters(ny+1:n, 2) = 0.0;   
                %   The exits are similar.
                exits(1:nx, 1) = uniform(obj, 0, 1, nx, 0);
                exits(1:nx, 2) = 1;
                exits(nx+1:n, 1) = 1;
                exits(nx+1:n, 2) = uniform(obj, 0, 1, ny, 1);
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