%> @file  Connect.m
%> @brief ERME connectivity matrix class definition.
% ==============================================================================
%> @brief Build a connectivity matrix for a Cartesian system.
% ==============================================================================
classdef Connect < handle
    
    properties (Access = private)
        %> User input.
        d_input
        %> Element map.
        d_elements
        d_number_elements
        d_M
        d_leak
    end
    
    methods 
        
        % ======================================================================
        %> @brief  Class constructor
        %
        %> Is is *assumed* that the element map is corrected.  The .m file based
        %> input has defined maps in real orientation; this has to be corrected
        %> during parsing.
        %>
        %> @param  input        Input database.
        %> @param  elements     Element map.
        %> @return              Instance of the Connect class.
        % ======================================================================
        function this = Connect(input, elements)
            this.d_input = input;
            this.d_elements = elements;
            this.d_number_elements = sum(sum(sum(elements>0)));
        end
        
        % ======================================================================
        %> @brief  Build a connectivity matrix.
        % ======================================================================
        function [M leak] = build(this)
            d = get(this.d_input, 'dimension');
            switch d
                case 1
                    build_1D(this);
                case 2
                    build_2D(this);
                case 3
                    build_3D(this);
                otherwise
                    error('Invalid dimension.');
            end
            M = this.d_M;
            leak = this.d_leak;
        end
        
        function plot_connect(this)
            spy(this.d_M>0, 'k');
            hold on;
            spy(this.d_M<0, 'r');
            hold off;
            title('connectivity matrix')
            xlabel(['nz = ',num2str(sum(sum(sum(this.d_M~=0))))]);
        end
        
    end
    
    methods  (Access = private)
        
        % ======================================================================
        %> @brief  Build a 1-d connectivity matrix.
        % ======================================================================
        function build_1D(this)
            
            bc_left  = get(this.d_input, 'bc_left');
            bc_right = get(this.d_input, 'bc_right');
            
            if sum(this.d_elements==0) > 0
                error('There are 0 elements in a slab problem!')
            end

            %-------------------------------------------------------------------
            % SETUP FOR CONNECTING FACES.  Matrix "mm" connects faces to faces 
            %   (or to global boundaries) with
            %      value of 1 = reflective
            %      value of 2 = connection between faces of different elements
            %   Note: all the "2"s due to 2 faces per element
            mm = zeros( 2*this.d_number_elements );
            I = this.d_number_elements;
            
            for i = 1:I
                % Connect my right to neighbor's left
                if i < I && this.d_elements(i+1) > 0
                    mm(2*(i-1)+2, 2*i+1) = 2;
                end
                % same for someone to my left
                if i > 1 && this.d_elements(i-1) > 0
                    mm(2*(i-1)+1, 2*(i-2)+2 ) = 2;
                end
                
                % SIMPLE TREATMENT OF BOUNDARIES
                if i == 1 % LEFT BOUNDARY
                    if strcmp(bc_left, 'reflect')
                        mm(2*(i-1)+1, 2*(i-1)+1 ) = 1;
                    end % else it stays zero
                end
                if i == I % RIGHT BOUNDARY
                    if strcmp(bc_right, 'reflect')
                        mm(2*(i-1)+2, 2*(i-1)+2 ) = 1;
                    end % else it stays zero
                end
            end
            mm = sparse(mm);
            
            %-------------------------------------------------------------------
            % SETUP FOR REFLECTION.  Because current responses are computed 
            %  w/r to the INCOMING coordinates, the outgoing response is 
            %  exactly the input of the neighbor cell; however, because 
            %  reflection switches left-to-right to right-to-left, all the 
            %  odd order legendre moments must have their
            %  signs reversed. That's why there is the -1 in the definition.
            ng = get(this.d_input, 'number_groups');
            po = get(this.d_input, 'rf_order_polar');  % space/azimuth don't matter    
            tot = (po+1);
            refl_v = zeros(tot, 1);
            i = 0;
            
            for p = 0:po
                i = i + 1;
                refl_v(i) = (-1)^mod(p, 2);
            end
            
            order = tot;
            refl_v_len = length(refl_v);
            refl  = spdiags( refl_v, 0, refl_v_len, refl_v_len);
            refl  = kron(  speye(ng), refl );
            
            %-------------------------------------------------------------------
            % SETUP FOR TRANSMISSION.  Because the current responses are 
            %  defined as input for adjacent cells, we need only a 1-to-1 map, 
            %  and so for each fact-to-face location (denoted by a 2 in mm), 
            %  we place a ones vector of appropriate length (i.e. # orders x 
            %  # groups), but in matrix form, i.e. an identity matrix
            tran = speye(order);
            tran = kron(  speye(ng), tran );
            
            %-------------------------------------------------------------------
            % COMBINE EVERYTHING IN M
            M   = kron(mm==1, refl )+kron(mm==2, tran);
            
        end
        
        % ======================================================================
        %> @brief  Build a 2-d connectivity matrix.
        % ======================================================================
        function build_2D(this)
            
            bc_left   = get(this.d_input, 'bc_left');
            bc_right  = get(this.d_input, 'bc_right');
            bc_bottom = get(this.d_input, 'bc_bottom');
            bc_top    = get(this.d_input, 'bc_top');

            %-------------------------------------------------------------------
            % SETUP FOR CONNECTING FACES.  Matrix "mm" connects faces to faces 
            %   (or to global boundaries) with
            %      value of 1 = reflective
            %      value of 2 = connection between faces of different elements
            %   Note: all the "4"s due to four faces per element
            mm = zeros( 4*this.d_number_elements );
            I = length(this.d_elements(:, 1));
            J = length(this.d_elements(1, :));
            k = 0;
            
            for j = 1:J      % rows
                for i = 1:I  % cols
                    if this.d_elements(i, j) > 0
                        % I am the kth element
                        k = k + 1;
                    else
                        % else I'm just void
                        continue
                    end
                    % if there is somebody (>0) to my right
                    % connect my face 2 to their face 1
                    if i < I && this.d_elements(i+1,j) > 0
                        mm( 4*(k-1)+2,4*k+1 ) = 2;
                    end
                    % same for someone to my left
                    if i > 1 && this.d_elements(i-1,j) > 0
                        mm( 4*(k-1)+1,4*(k-2)+2 ) = 2;
                    end
                    % if there is somebody (>0) to my top
                    % connect my face 4 to their face 3
                    if j < J && this.d_elements(i,j+1) > 0
                        % nzrow counts # elements left in my row and in
                        % next row before my neighbor
                        nzrow = sum( this.d_elements(i+1:end,j)>0 );
                        nzrow = nzrow + sum( this.d_elements(1:i,j+1)>0 );
                        mm( 4*(k-1)+4, 4*(k-1+nzrow)+3 ) = 2;
                    end
                    % same for someone to my bottom
                    if j > 1 && this.d_elements(i,j-1) > 0
                        nzrow = sum( this.d_elements(1:i,j)>0 );
                        nzrow = nzrow + sum( this.d_elements(i+1:end,j-1)>0 );
                        mm( 4*(k-1)+3,4*(k-1-nzrow)+4 ) = 2;
                    end
                    % SIMPLE TREATMENT OF BOUNDARIES
                    if i == 1 % LEFT BOUNDARY
                        if strcmp(bc_left, 'reflect')
                            mm( 4*(k-1)+1, 4*(k-1)+1 ) = 1;
                        else
                            mm( 4*(k-1)+1, 4*(k-1)+1 ) = -1;
                        end
                    end
                    if i == I % RIGHT BOUNDARY
                        if strcmp(bc_right, 'reflect')
                            mm( 4*(k-1)+2, 4*(k-1)+2 ) = 1;
                        else
                            mm( 4*(k-1)+2, 4*(k-1)+2 ) = -1;
                        end
                    end
                    if j == 1 % BOTTOM BOUNDARY
                        if strcmp(bc_bottom, 'reflect')
                            mm( 4*(k-1)+3, 4*(k-1)+3 ) = 1;
                        else
                            mm( 4*(k-1)+3, 4*(k-1)+3 ) = -1;
                        end
                    end
                    if j == J % TOP BOUNDARY
                        if strcmp(bc_top, 'reflect')
                            mm( 4*(k-1)+4, 4*(k-1)+4 ) = 1;
                        else
                            mm( 4*(k-1)+4, 4*(k-1)+4 ) = -1;
                        end
                    end
                end
            end
            mm = sparse(mm);
            
            
            %-------------------------------------------------------------------
            % SETUP FOR REFLECTION.  Because current responses are computed 
            %  w/r to the INCOMING coordinates, the outgoing response is 
            %  exactly the input of the neighbor cell; however, because 
            %  reflection switches left-to-right to right-to-left, all the 
            %  odd order legendre moments must have their
            %  signs reversed. That's why there is the -1 in the definition.
            ng = get(this.d_input, 'number_groups');
            so = get(this.d_input, 'rf_order_space');
            ao = get(this.d_input, 'rf_order_azimuth');
            po = get(this.d_input, 'rf_order_polar');       
            tot = (so+1)*(ao+1)*(po+1);
            refl_v = zeros(tot, 1);
            i = 0;
            for s = 0:so
                for a = 0:ao
                    for p = 0:po
                        i = i + 1;
                        refl_v(i) = (-1)^mod(s+a, 2);
                    end
                end
            end
            order = tot;
            refl_v_len = length(refl_v);
            refl  = spdiags( refl_v, 0, refl_v_len, refl_v_len);
            refl  = kron(  speye(ng), refl );
            
            %-------------------------------------------------------------------
            % SETUP FOR TRANSMISSION.  Because the current responses are 
            %  defined as input for adjacent cells, we need only a 1-to-1 map, 
            %  and so for each fact-to-face location (denoted by a 2 in mm), 
            %  we place a ones vector of appropriate length (i.e. # orders x 
            %  # groups), but in matrix form, i.e. an identity matrix
            tran = speye(order);
            tran = kron(  speye(ng), tran );
            
            %-------------------------------------------------------------------
            % COMBINE EVERYTHING IN M
            this.d_M = kron(mm==1, refl )+kron(mm==2, tran);
            this.d_leak = ((mm==-1)*ones(4*this.d_number_elements, 1))';
            % Leakage vector used as follows:
            %  leakage = leak * (L*J)
            
        end
        
        % ======================================================================
        %> @brief  Build a 3-d connectivity matrix.
        % ======================================================================
        function build_3D(this)
            
            bc_left   = get(this.d_input, 'bc_left');
            bc_right  = get(this.d_input, 'bc_right');
            bc_bottom = get(this.d_input, 'bc_bottom');
            bc_top    = get(this.d_input, 'bc_top');
            bc_north  = get(this.d_input, 'bc_north');
            bc_south  = get(this.d_input, 'bc_south');

            %-------------------------------------------------------------------
            % SETUP FOR CONNECTING FACES.  Matrix "mm" connects faces to faces 
            %   (or to global boundaries) with
            %      value of 1 = reflective
            %      value of 2 = connection between faces of different elements
            %   Note: all the "6"s due to six faces per element
            mm = zeros( 6*this.d_number_elements );
            I = length(this.d_elements(:, 1, 1));
            J = length(this.d_elements(1, :, 1));
            K = length(this.d_elements(1, 1, :));
            
            e = 0
            for k = 1:K
            for j = 1:J      
                for i = 1:I  
                    if this.d_elements(i, j, k) > 0
                        % I am the eth element
                        e = e + 1;
                    else
                        % else I'm just void
                        continue
                    end
                    
                    % if there is somebody (>0) to my right
                    % connect my face 2 to their face 1
                    if i < I && this.d_elements(i+1, j, k) > 0
                        mm( 6*(e-1)+2, 6*e+1 ) = 1;
                    end
                    % same for someone to my left
                    if i > 1 && this.d_elements(i-1, j, k) > 0
                        mm( 6*(e-1)+1, 6*(e-2)+2 ) = 2;
                    end
                    
                    % if there is somebody (>0) to my top
                    % connect my face 4 to their face 3
                    if j < J && this.d_elements(i, j+1, k) > 0
                        % nzrow counts # elements left in my row and in
                        %   next row before my neighbor
                        nzrow = sum( this.d_elements(i+1:end, j, k)>0 );
                        nzrow = nzrow + sum( this.d_elements(1:i, j+1, k)>0 );
                        mm( 6*(e-1)+4, 6*(e-1+nzrow)+3 ) = 3;
                    end
                    % same for someone to my bottom
                    if j > 1 && this.d_elements(i,j-1,k) > 0
                        nzrow = sum( this.d_elements(1:i, j, k)>0 );
                        nzrow = nzrow + sum( this.d_elements(i+1:end, j-1, k)>0 );
                        mm( 6*(e-1)+3, 6*(e-1-nzrow)+4 ) = 4;
                    end
                    
                    % if there is somebody (>0) to my north
                    % connect my face 6 to their face 5
                    if k < K && this.d_elements(i, j, k+1) > 0
                        % nzrow counts # elements left in my row and in
                        %   next row before my neighbor
                        nzrow = sum( this.d_elements(i+1:end, j, k)>0 );              % rest if this row
                        nzrow = nzrow + sum(sum( this.d_elements(:, j+1:end, k)>0 )); % rest of this plane
                        nzrow = nzrow + sum(sum( this.d_elements(:, 1:j-1, k+1)>0 )); % plane before this row one level up
                        nzrow = nzrow + sum( this.d_elements(1:i, j, k+1)>0 );        % up to this point one level up
                        mm( 6*(e-1)+6, 6*(e-1+nzrow)+5 ) = 5;
                    end
                    % same for someone to my bottom
                    if k > 1 && this.d_elements(i, j, k-1) > 0
                        nzrow = sum( this.d_elements(1:i, j, k)>0 );
                        nzrow = nzrow + sum(sum( this.d_elements(:, 1:j-1, k)>0 ));
                        nzrow = nzrow + sum(sum( this.d_elements(:, j+1:end, k-1)>0 ));
                        nzrow = nzrow + sum( this.d_elements(i+1:end, j, k-1)>0 );
                        mm( 6*(e-1)+5, 6*(e-1-nzrow)+6 ) = 6;
                    end                    
                    
                    % SIMPLE TREATMENT OF BOUNDARIES
                    if i == 1 % LEFT BOUNDARY
                        if strcmp(bc_left, 'reflect')
                            mm( 6*(e-1)+1, 6*(e-1)+1 ) = -1;
                        end % else it stays zero
                    end
                    if i == I % RIGHT BOUNDARY
                        if strcmp(bc_right, 'reflect')
                            mm( 6*(e-1)+2, 6*(e-1)+2 ) = -1;
                        end % else it stays zero
                    end
                    if j == 1 % BOTTOM BOUNDARY
                        if strcmp(bc_bottom, 'reflect')
                            mm( 6*(e-1)+3, 6*(e-1)+3 ) = -1;
                        end % else it stays zero
                    end
                    if j == J % TOP BOUNDARY
                        if strcmp(bc_top, 'reflect')
                            mm( 6*(e-1)+4, 6*(e-1)+4 ) = -1;
                        end % else it stays zeros
                    end
                    if k == 1 % SOUTH BOUNDARY
                        if strcmp(bc_south, 'reflect')
                            mm( 6*(e-1)+5, 6*(e-1)+5 ) = -2;
                        end % else it stays zero
                    end
                    if k == K % NORTH BOUNDARY
                        if strcmp(bc_north, 'reflect')
                            mm( 6*(e-1)+6, 6*(e-1)+6 ) = -2;
                        end % else it stays zeros
                    end                    
                end
            end
            end
            mm = sparse(mm);
            
            
            %-------------------------------------------------------------------
            % SETUP FOR REFLECTION.  Because current responses are computed 
            %  w/r to the INCOMING coordinates, the outgoing response is 
            %  exactly the input of the neighbor cell; however, because 
            %  reflection switches left-to-right to right-to-left, all the 
            %  odd order legendre moments must have their
            %  signs reversed. That's why there is the -1 in the definition.
            ng = get(this.d_input, 'number_groups');
            so = get(this.d_input, 'rf_order_space');
            ao = get(this.d_input, 'rf_order_azimuth');
            po = get(this.d_input, 'rf_order_polar');       
            tot = (so+1)*(so+1)*(ao+1)*(po+1);
            refl_V = zeros(tot, 1); % For the vertical faces
            refl_H = zeros(tot, 1); % For south and north faces (horizontal)
            i = 0;
            for sy = 0:so
                for sx = 0:so
                    for a = 0:ao
                        for p = 0:po
                            i = i + 1;
                            refl_V(i) = (-1)^mod(sx+sy+a,   2); % A particle going up stays going up off a V face-->no polar polarity
                            refl_H(i) = (-1)^mod(sx+sy+a+p, 2); % Off north/south, there *is* polar polarity switch
                        end
                    end
                end
            end
            order = tot;
            refl_len = length(refl_V);
            refl_V  = spdiags( refl_V, 0, refl_len, refl_len);
            refl_V  = kron(  speye(ng), refl_V );
            refl_H  = spdiags( refl_H, 0, refl_len, refl_len);
            refl_H  = kron(  speye(ng), refl_H );            
            
            %-------------------------------------------------------------------
            % SETUP FOR TRANSMISSION.  Because the current responses are 
            %  defined as input for adjacent cells, we need only a 1-to-1 map, 
            %  and so for each fact-to-face location (denoted by a 2 in mm), 
            %  we place a ones vector of appropriate length (i.e. # orders x 
            %  # groups), but in matrix form, i.e. an identity matrix
            tran = speye(order);
            tran = kron(  speye(ng), tran );
            
            %-------------------------------------------------------------------
            % COMBINE EVERYTHING IN M
            this.d_M = kron(mm==-2, refl_H )+kron(mm==-1, refl_V )+kron(mm>0, tran); 
        end
        
    end
    
end




