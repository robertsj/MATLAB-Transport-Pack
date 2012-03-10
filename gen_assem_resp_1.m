%> @file  test_response.m
%> @brief Tests a response function boundary problem.
%>
%> The problem is a two group node with multiplication (that can be turned
%> on or off).  We're looking for the outgoing fluxes.
%>
%> Issues:
%>     1) Look at using Krylov for the entire group solve
%>     2) GMRES for response (bc's not updated right now...)
%>     3) Multigroup acceleration?  CMFD?
%>     4) RF table format 1.343001
% there is
%
% Add the capability to compute the flux for a given set of incident responses.
% That way, we can see what higher current modes might mean physically.
%
% Also, does the "adjoint" fall out somewhere?
% ==============================================================================
%clear
flag = 1;

%clear classes

% Get the default input.
input = Input();

% Material etc.
ng = 2;
put(input, 'number_groups',         ng);

% Inner iteration parameters.
put(input, 'eigen_tolerance',       1e-9);
put(input, 'eigen_max_iters',       2000);

put(input, 'inner_tolerance',       1e-9);
put(input, 'inner_max_iters',       100);
put(input, 'outer_tolerance',       1e-9);
put(input, 'outer_max_iters',       10);
put(input, 'inner_solver',          'SI');
put(input, 'livolant_free_iters',   3);
put(input, 'livolant_accel_iters',  6);
if flag == 0
put(input, 'bc_left',               'reflect');
put(input, 'bc_right',              'reflect');
put(input, 'bc_top',                'reflect');
put(input, 'bc_bottom',             'reflect');
else
put(input, 'bc_left',               'response');
put(input, 'bc_right',              'vacuum');
put(input, 'bc_top',                'vacuum');
put(input, 'bc_bottom',             'vacuum');
end
put(input, 'print_out',             1);

% Set the incident response order
put(input, 'rf_order_group',        0);
put(input, 'rf_order_space',        0);
put(input, 'rf_order_polar',        0);
put(input, 'rf_order_azimuth',      0);

put(input, 'rf_max_order_space',        2);
put(input, 'rf_max_order_azimuth',      2);
put(input, 'rf_max_order_polar',        0);

qr=2;
quadrature  = QuadrupleRange(qr);
if qr==2
put(input, 'quad_number_polar',     1);
put(input, 'quad_number_azimuth',   2);
elseif qr==8
put(input, 'quad_number_polar',     2);
put(input, 'quad_number_azimuth',   4);
elseif qr==18
put(input, 'quad_number_polar',     3);
put(input, 'quad_number_azimuth',   6);
else
   error('angle wrong') 
end


elements = [1];
number_elements = 1;        

% elements = [1 ];
% number_elements = 1; 

M = Connect(input, elements, number_elements);


% ==============================================================================
% MATERIALS (Test two group data)
% ==============================================================================
mat = test_materials(2);

tic
% ==============================================================================
% PINS
% ==============================================================================

% Shared pin properties
pitch = 1.26; radii = 0.54; number = 10;
% Pin 1 - Fuel 1
matid = [2 1];  pin1 = PinCell(pitch, radii, matid); meshify(pin1, number);
% Pin 2 - Fuel 2
matid = [2 1];  pin2 = PinCell(pitch, radii, matid); meshify(pin2, number);
% Pin 3 - MOD
matid = [  2];  pin3 = PinCell(pitch,    [], matid); meshify(pin3, number);

% ==============================================================================
% ASSEMBLIES 0.701423630411957
% ==============================================================================

% Assembly 1 
pin_map1 = [1 1 1
            1 2 1
            1 1 1];     
mesh = Assembly({pin1, pin2, pin3}, pin_map1); 
meshify(mesh);
toc

% ==============================================================================
% SETUP
% ==============================================================================
tic
state       = State(input, mesh);
boundary    = BoundaryMesh(input, mesh, quadrature);
% Empty external source.
q_e         = Source(mesh, ng);
% Use a fission source.
q_f = FissionSource(state, mesh, mat);
initialize(q_f);
toc
if (flag==0)
% Make the inner iteration.
solver= Eigensolver(input,        ...
              state,    	...
              boundary,     ...
              mesh,     	...
              mat,        	...
              quadrature, 	...
              q_e,          ...
              q_f);
 solve(solver);          %  1.27828651128439
error('stop')
end

% RESPONSE LOOP
k = 0;
coef = cell(4, 1);
max_s_o = get(input, 'rf_max_order_space'); 
max_a_o = get(input, 'rf_max_order_azimuth');
max_p_o = get(input, 'rf_max_order_polar');


max_o   = max_s_o * max_a_o * max_p_o;
coef{1} = zeros(max_o);
coef{2} = zeros(max_o);
coef{3} = zeros(max_o);
coef{4} = zeros(max_o);
total = ng*(1+max_s_o)*(1+max_a_o)*(1+max_p_o);

% Get the left bc
boundary = BoundaryMesh(input, mesh, quadrature);
bc = get_bc(boundary, Mesh.LEFT);

solver= KrylovMG(   input,        ...
                    state,    	...
                    boundary,     ...
                    mesh,     	...
                    mat,        	...
                    quadrature, 	...
                    q_e,          ...
                    q_f);


t0 = tic;
t1 = t0;
t2 = 0;
for g_o = 1:ng
for s_o = 0:max_s_o
    for a_o = 0:max_a_o
        for p_o = 0:max_p_o
            
            k = k + 1;  
            disp([' doing ',num2str(k),' of ',num2str(total)])
            if (k > 1)
                elap = toc;
                per  = elap / k;
                tot  = per * (total-k);
                disp([' estimated remaining time = ',num2str(tot),' seconds.']) 
            end
            set_orders(bc, g_o, s_o, a_o, p_o);
            
            set_keff(solver,  1.325242580556770); % kinf 0.9009968666

            out = solve(solver);
            reset(q_f);

            % make basis
            number_space   = number_cells_x(mesh);
            number_polar   = get(input, 'quad_number_polar');
            number_azimuth = get(input, 'quad_number_azimuth');                
            basis_space    = DiscreteLP(number_space-1);
            basis_polar    = DiscreteLP(number_polar-1);
            basis_azimuth  = DiscreteLP(number_azimuth*2-1);
            num_ang        = number_polar*number_azimuth;

            octants = [3 2   % ref
                       1 4   % far
                       4 3   % lef
                       2 1]; % rig

            % Expand the coefficients
            for side = 1:4
                            
                % always left to right in space w/r to outgoing flux
                if side > 0 || side == 2 || side == 3
                    s1 = 1;
                    s2 = number_space;
                    s3 = 1;
                else
                    s1 = number_space;
                    s2 = 1;
                    s3 = -1;
                end
                
                for gg = 1:ng                    
                    for o = 1:2
                        o_in  = octants(side, o); % incident octant
                        % always left to right in angle w/r to outgoing flux
                        if o == 1
                            a1 = 1;
                            a2 = number_azimuth;
                            a3 = 1;
                        else
                            a1 = number_azimuth;
                            a2 = 1;
                            a3 = -1;
                        end
                        % Get psi(x, angles)
                        initialize(boundary, gg);
                        if (side == 1 || side == 2)
                            ft = get_psi_v_octant(boundary, o_in, Boundary.OUT);
                            if (side == 1)
                                ft(1:end,:)=ft(end:-1:1,:); % reverse space
                            end
                        else
                            ft = get_psi_h_octant(boundary, o_in, Boundary.OUT);
                            if (side == 4)
                                ft(1:end,:)=ft(end:-1:1,:);  % reverse space
                            end
                        end
                        % populate the vectors we expand
                        for s = s1:s3:s2
                            ang = 1;
                            for a = a1:a3:a2
                                for p = 1:number_polar
                                    f(s, (o-1)*num_ang+ang,gg) = ft(s, (a-1)*number_polar+p);
                                    ang = ang + 1;
                                end
                            end
                        end
                    end
                    
                end

                % Group->Space->Azimuth->Polar.  f(space, angle, group)
                i = 1;
                for gg = 1:ng
                for ord_s = 1:max_s_o+1
                    for ord_a = 1:max_a_o+1
                        for ord_p = 1:max_p_o+1
                            psi_ap = zeros(number_azimuth, number_polar);
                            angle = 0;
                            az=0;
                            for o = 1:2
                            if (side == Mesh.LEFT || side == Mesh.RIGHT)
                                oo = 1;
                            else
                                oo = 2;
                            end
                            if (1 == 1) % this makes angles symmetric
                                a1 = 1; 
                                a2 = number_azimuth;
                                a3 = 1;
                            else
                                a1 = number_azimuth;
                                a2 = 1;
                                a3 = -1;
                            end
                            
                            for a = a1:a3:a2
                                az = az+1;
                                for p = 1:number_polar
                                    angle = angle + 1;
                                    psi_ap(az, p) = f(:, angle,gg)'*basis_space(:, ord_s);
                                end
                            end
                            end
                            
                            psi_p = zeros(number_polar, 1); 
                            b = basis_azimuth(:, ord_a);
                            if side > 2
                                lb = length(b);
                                b(1:lb/2) = b(lb/2:-1:1);
                                b(lb/2+1:end) = b(end:-1:lb/2+1);
                            end
                            for p = 1:number_polar
                                psi_p(p) = psi_ap(:, p)'*b;
                            end
                            coef{side}(i, k) = psi_p'*basis_polar(:, ord_p); % i <- k
                            i = i + 1;
                        end
                    end
                end
                end
                
            end % side loop
            
        end % azimuth loop
    end % polar loop
end % space loop
end % group loop
toc
clear R
max_o = ng*(max_s_o + 1)*(max_a_o + 1)*(max_p_o + 1);
R( (0*max_o)+1: 1*max_o, (0*max_o)+1: 1*max_o) = coef{1}(:, :); % left -> left
R( (1*max_o)+1: 2*max_o, (0*max_o)+1: 1*max_o) = coef{2}(:, :); % left -> right
R( (2*max_o)+1: 3*max_o, (0*max_o)+1: 1*max_o) = coef{3}(:, :); % left -> top
R( (3*max_o)+1: 4*max_o, (0*max_o)+1: 1*max_o) = coef{4}(:, :); % left -> bottom


R( (0*max_o)+1: 1*max_o, (1*max_o)+1: 2*max_o) = coef{2}(:, :); % right -> lef
R( (1*max_o)+1: 2*max_o, (1*max_o)+1: 2*max_o) = coef{1}(:, :); % right -> right
R( (2*max_o)+1: 3*max_o, (1*max_o)+1: 2*max_o) = coef{4}(:, :); % right -> top
R( (3*max_o)+1: 4*max_o, (1*max_o)+1: 2*max_o) = coef{3}(:, :); % right -> bottom


R( (0*max_o)+1: 1*max_o, (2*max_o)+1: 3*max_o) = coef{4}(:, :); % bottom -> left
R( (1*max_o)+1: 2*max_o, (2*max_o)+1: 3*max_o) = coef{3}(:, :); % bottom -> right
R( (2*max_o)+1: 3*max_o, (2*max_o)+1: 3*max_o) = coef{1}(:, :); % bottom -> bottom
R( (3*max_o)+1: 4*max_o, (2*max_o)+1: 3*max_o) = coef{2}(:, :); % bottom -> top


R( (0*max_o)+1: 1*max_o, (3*max_o)+1: 4*max_o) = coef{3}(:, :); % top -> left
R( (1*max_o)+1: 2*max_o, (3*max_o)+1: 4*max_o) = coef{4}(:, :); % top -> right
R( (2*max_o)+1: 3*max_o, (3*max_o)+1: 4*max_o) = coef{2}(:, :); % top -> bottom
R( (3*max_o)+1: 4*max_o, (3*max_o)+1: 4*max_o) = coef{1}(:, :); % top -> top


RR = kron(speye(number_elements), R);

e = eigs(M*RR, 4, 'LR')

lambda(MSO+1,MAO+1,MPO+1) = e(1);

