function M = Connect(input, elements, number_elements)
% function M = Connect(input, elements, number_elements)
%            input = matrix of element placements
%         elements = number of elements
%  number_elements = number of groups
%  
%  Basic extension of serment-m version.  Here, we handle higher order angular
%  expansions.

% Assign the boundary conditions
lef = 1;
rig = 0;
bot = 1;
top = 0;

%--------------------------------------------------------------------------
% SETUP FOR CONNECTING FACES.  Matrix "mm" connects faces to faces (or 
%   to global boundaries) with
%      value of 1 = reflective
%      value of 2 = connection between faces of different elements
%   Note: all the "4"s due to four faces per element
mm = zeros( 4*number_elements );
I = length(elements(:,1));
J = length(elements(1,:));
k = 0;

for j = 1:J      % rows
    for i = 1:I  % cols
        if elements(i,j) > 0
            % I am the kth element
            k = k + 1;
        else
            % else I'm just void
            continue
        end        
        % if there is somebody (>0) to my right
        % connect my face 2 to their face 1
        if i < I && elements(i+1,j) > 0 
            mm( 4*(k-1)+2,4*k+1 ) = 2;
        end
        % same for someone to my left
        if i > 1 && elements(i-1,j) > 0
            mm( 4*(k-1)+1,4*(k-2)+2 ) = 2;
        end            
        % if there is somebody (>0) to my top
        % connect my face 4 to their face 3
        if j < J && elements(i,j+1) > 0 
            % nzrow counts # elements left in my row and in
            % next row before my neighbor
            nzrow = sum( elements(i+1:end,j)>0 );
            nzrow = nzrow + sum( elements(1:i,j+1)>0 );
            mm( 4*(k-1)+4, 4*(k-1+nzrow)+3 ) = 2;
        end
        % same for someone to my bottom
        if j > 1 && elements(i,j-1) > 0
            nzrow = sum( elements(1:i,j)>0 );
            nzrow = nzrow + sum( elements(i+1:end,j-1)>0 );
            mm( 4*(k-1)+3,4*(k-1-nzrow)+4 ) = 2;
        end       
        % SIMPLE TREATMENT OF BOUNDARIES
        if i == 1 % LEFT BOUNDARY
            if lef == 1
                mm( 4*(k-1)+1, 4*(k-1)+1 ) = 1;
            end % else it stays zero
        end
        if i == I % RIGHT BOUNDARY
            if rig == 1
                mm( 4*(k-1)+2, 4*(k-1)+2 ) = 1;
            end % else it stays zero
        end
        if j == 1 % BOTTOM BOUNDARY
            if bot == 1
                mm( 4*(k-1)+3, 4*(k-1)+3 ) = 1;
            end % else it stays zero
        end    
        if j == J % TOP BOUNDARY
            if top == 1
                mm( 4*(k-1)+4, 4*(k-1)+4 ) = 1;
            end % else it stays zeros
        end            
    end
end
mm = sparse(mm);


%--------------------------------------------------------------------------
% SETUP FOR REFLECTION.  Because current responses are computed w/r to
%  the INCOMING coordinates, the outgoing response is exactly the input
%  of the neighbor cell; however, because reflection switches left-to-right
%  to right-to-left, all the odd order legendre moments must have their
%  signs reversed. That's why there is the -1 in the definition.
ng = get(input, 'number_groups');
so = get(input, 'rf_max_order_space');
ao = get(input, 'rf_max_order_azimuth');
po = get(input, 'rf_max_order_polar');
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

%refl_v     = ones(order+1,1).*(-1).^(1+(1:order+1)');
refl_v_len = length(refl_v);
refl  = spdiags( refl_v, 0, refl_v_len, refl_v_len);
refl  = kron(  speye(ng), refl );

%--------------------------------------------------------------------------
% SETUP FOR TRANSMISSION.  Because the current responses are defined as
%  input for adjacent cells, we need only a 1-to-1 map, and so for each
%  fact-to-face location (denoted by a 2 in mm), we place a ones vector
%  of appropriate length (i.e. # orders x # groups), but in matrix form,
%  i.e. an identity matrix
tran = speye(order);
tran = kron(  speye(ng), tran );

%--------------------------------------------------------------------------
% COMBINE EVERYTHING IN M
M   = kron(mm==1, refl )+kron(mm==2, tran);


