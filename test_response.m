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
%>     4) RF table format
% there is
%
% Add the capability to compute the flux for a given set of incident responses.
% That way, we can see what higher current modes might mean physically.
%
% Also, does the "adjoint" fall out somewhere?
% ==============================================================================

flag = 0;

%clear classes

% Get the default input.
input = Input();

% Material etc.
put(input, 'number_groups',         1);

% Inner iteration parameters.
put(input, 'inner_tolerance',       1e-10);
put(input, 'inner_max_iters',       100);
put(input, 'outer_tolerance',       1e-8);
put(input, 'outer_max_iters',       300);
put(input, 'inner_solver',          'SI');
put(input, 'livolant_free_iters',   3);
put(input, 'livolant_accel_iters',  6);
if flag == 0
put(input, 'bc_left',               'reflect');
put(input, 'bc_right',              'vacuum');
put(input, 'bc_top',                'vacuum');
put(input, 'bc_bottom',             'reflect');
else
put(input, 'bc_left',               'response');
put(input, 'bc_right',              'vacuum');
put(input, 'bc_top',                'vacuum');
put(input, 'bc_bottom',             'vacuum');
end
put(input, 'print_out',             1);

% Set the incident response order
put(input, 'rf_order_group',        1);
put(input, 'rf_order_space',        0);
put(input, 'rf_order_polar',        0);
put(input, 'rf_order_azimuth',      0);
put(input, 'rf_max_order_space',        2);
put(input, 'rf_max_order_azimuth',      3);
put(input, 'rf_max_order_polar',        0);

put(input, 'quad_number_polar',     1);
put(input, 'quad_number_azimuth',   2);


elements = [1 
            ];
number_elements = 1;        

M = Connect(input, elements, number_elements);


% One material, one group, c = 0.9
mat         = test_materials(1);

% Simple uniformly meshed square
mesh        = test_mesh(1);

state       = State(input, mesh);
quadrature  = QuadrupleRange(2);

boundary    = BoundaryMesh(input, mesh, quadrature);

% Empty external source.
q_e         = Source(mesh, 2);

% Use a fission source.
q_f = FissionSource(state, mesh, mat);
initialize(q_f);

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
 
%set_keff(solver, 0.49);
solve(solver);          %0.421086  0.5=0.99997268
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
total = (1+max_s_o)*(1+max_a_o)*(1+max_p_o);

for s_o = 0:max_s_o
    for a_o = 0:max_a_o
        for p_o = 0:max_p_o
            
            k = k + 1;  
            disp([' doing ',num2str(k),' of ',num2str(total)])
            % Set the incident response order
            put(input, 'rf_order_group',        1);
            put(input, 'rf_order_space',        s_o);
            put(input, 'rf_order_polar',        p_o);
            put(input, 'rf_order_azimuth',      a_o);

            boundary    = BoundaryMesh(input, mesh, quadrature);
            
            % Make the inner iteration. KrylovMG  FixedMultiply
            solver= KrylovMG(input,        ...
                state,    	...
                boundary,     ...
                mesh,     	...
                mat,        	...
                quadrature, 	...
                q_e,          ...
                q_f);
            set_keff(solver, 0.488); %491
            %solver.d_keff = 0.5;
            % Solve the problem
            tic
                out = solve(solver); 
            toc
            reset(q_f);
            %Get the flux
            phi{k} = flux(state, 1);
            %subplot(2,1,1)
           % figure(1)
           % plot_flux(mesh, phi{k})
           %error('fuck')
            axis square
            shading flat
            % Go through and get boundary responses

            % We're incident on the bottom.  
            %  ref -- bottom
            %  far -- top
            %  lef -- left
            %  rig -- right

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
                for o = 1:2
                    o_in  = octants(side, o); % incident octant
                    % always left to right in angle w/r to outgoing flux
                    if o == 1
                        a1 = 1;
                        a2 = num_ang;
                        a3 = 1;
                    else
                        a1 = num_ang;
                        a2 = 1;
                        a3 = -1;
                    end
                    % Get psi(x, angles)
                    if (side == 1 || side == 2)
                        ft = get_psi_v_octant(boundary, o_in, Boundary.OUT);
                        if (side == 1)
                           ft(1:end,:)=ft(end:-1:1,:); 
                        end
                    else
                        ft = get_psi_h_octant(boundary, o_in, Boundary.OUT);  
                        if (side == 4)
                           ft(1:end,:)=ft(end:-1:1,:); 
                        end                        
                    end
                    for s = s1:s3:s2
                        ang = 1;
                        for a = a1:a3:a2
                            f(s, (o-1)*num_ang+ang) = ft(s, a);
                            ang = ang + 1;
                        end
                    end
                end

                % Space->Azimuth->Polar.  f(space, angle, group)
                i = 1;
                for ord_s = 1:max_s_o+1
                    for ord_a = 1:max_a_o+1
                        for ord_p = 1:max_p_o+1
                            psi_ap = zeros(number_azimuth, number_polar);
                            angle = 0;
                            for a = 1:number_azimuth*2
                                for p = 1:number_polar
                                    angle = angle + 1;
                                    psi_ap(a, p) = f(:, angle)'*basis_space(:, ord_s);
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
                
            end % side loop
            
        end % azimuth loop
    end % polar loop
end % space loop

max_o = (max_s_o + 1)*(max_a_o + 1)*(max_p_o + 1);
R( (0*max_o)+1: 1*max_o, (0*max_o)+1: 1*max_o) = coef{1}(:, :); % left -> left
R( (1*max_o)+1: 2*max_o, (0*max_o)+1: 1*max_o) = coef{2}(:, :); % left -> right
R( (2*max_o)+1: 3*max_o, (0*max_o)+1: 1*max_o) = coef{3}(:, :); % left -> top
R( (3*max_o)+1: 4*max_o, (0*max_o)+1: 1*max_o) = coef{4}(:, :); % left -> bottom


R( (0*max_o)+1: 1*max_o, (1*max_o)+1: 2*max_o) = coef{2}(:, :); % right -> lef
R( (1*max_o)+1: 2*max_o, (1*max_o)+1: 2*max_o) = coef{1}(:, :); % right -> right
R( (2*max_o)+1: 3*max_o, (1*max_o)+1: 2*max_o) = coef{4}(:, :); % right -> top
R( (3*max_o)+1: 4*max_o, (1*max_o)+1: 2*max_o) = coef{3}(:, :); % right -> bottom


R( (0*max_o)+1: 1*max_o, (2*max_o)+1: 3*max_o) = coef{3}(:, :); % bottom -> left
R( (1*max_o)+1: 2*max_o, (2*max_o)+1: 3*max_o) = coef{4}(:, :); % bottom -> right
R( (2*max_o)+1: 3*max_o, (2*max_o)+1: 3*max_o) = coef{1}(:, :); % bottom -> bottom
R( (3*max_o)+1: 4*max_o, (2*max_o)+1: 3*max_o) = coef{2}(:, :); % bottom -> top


R( (0*max_o)+1: 1*max_o, (3*max_o)+1: 4*max_o) = coef{4}(:, :); % top -> left
R( (1*max_o)+1: 2*max_o, (3*max_o)+1: 4*max_o) = coef{3}(:, :); % top -> right
R( (2*max_o)+1: 3*max_o, (3*max_o)+1: 4*max_o) = coef{2}(:, :); % top -> bottom
R( (3*max_o)+1: 4*max_o, (3*max_o)+1: 4*max_o) = coef{1}(:, :); % top -> top


RR = kron(speye(number_elements), R);

eigs(M*RR)
