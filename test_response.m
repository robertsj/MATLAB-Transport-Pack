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
% ==============================================================================



%clear classes

% Get the default input.
input = Input();

% Material etc.
put(input, 'number_groups',         1);

% Inner iteration parameters.
put(input, 'inner_tolerance',       1e-4);
put(input, 'inner_max_iters',       100);
put(input, 'outer_tolerance',       1e-4);
put(input, 'outer_max_iters',       200);
put(input, 'inner_solver',          'Livolant');
put(input, 'livolant_free_iters',   3);
put(input, 'livolant_accel_iters',  6);
put(input, 'bc_left',               'vacuum');
put(input, 'bc_right',              'vacuum');
put(input, 'bc_top',                'vacuum');
put(input, 'bc_bottom',             'response');
put(input, 'print_out',             1);

% Set the incident response order
put(input, 'rf_order_group',        1);
put(input, 'rf_order_space',        0);
put(input, 'rf_order_polar',        0);
put(input, 'rf_order_azimuth',      0);
put(input, 'quad_number_polar',     1);
put(input, 'quad_number_azimuth',   2);

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

% Make the inner iteration.
solver= FixedMultiply(input,        ...
              state,    	...
              boundary,     ...
              mesh,     	...
              mat,        	...
              quadrature, 	...
              q_e,          ...
              q_f);
 
          
          
% RESPONSE LOOP
k = 0;
coef = cell(4, 1);
coef{1} = zeros(2*2*1, 2*2*1);
coef{2} = zeros(2*2*1, 2*2*1);
coef{3} = zeros(2*2*1, 2*2*1);
coef{4} = zeros(2*2*1, 2*2*1);
for s_o = 0:1
    for a_o = 0:1
        for p_o = 0:0
            k = k + 1;       
            % Set the incident response order
            put(input, 'rf_order_group',        1);
            put(input, 'rf_order_space',        s_o);
            put(input, 'rf_order_polar',        p_o);
            put(input, 'rf_order_azimuth',      a_o);

            boundary    = BoundaryMesh(input, mesh, quadrature);
            
            % Make the inner iteration.
            solver= FixedMultiply(input,        ...
                state,    	...
                boundary,     ...
                mesh,     	...
                mat,        	...
                quadrature, 	...
                q_e,          ...
                q_f);
             
            % Solve the problem
            tic
                out = solve(solver); 
            toc

            %Get the flux
            phi = flux(state, 1);
            %subplot(2,1,1)
            plot_flux(mesh, phi)
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


            octants = [4 3   % ref
                       2 1   % far
                       3 2   % lef
                       1 4]; % rig

            

            % Expand the coefficients
            for side = 1:4
                % always left to right in space w/r to outgoing flux
                if side == 2 || side == 3
                    s1 = 1;
                    s2 = number_space;
                    s3 = 1;
                else
                    s1 = number_space;
                    s2 = -1;
                    s3 = 1;
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
                    ft = get_psi_h_octant(boundary, o_in, Boundary.OUT);
                    for s = s1:s3:s2
                        angle = 1;
                        for a = a1:a3:a2
                            f(s, (o-1)*num_ang+angle) = ft(s, a);
                            angle = angle + 1;
                        end
                    end
                end

                % Space->Azimuth->Polar.  f(space, angle, group)
                i = 1;
                for ord_s = 1:2
                    for ord_a = 1:2
                        for ord_p = 1:1
                            psi_ap = zeros(number_azimuth, number_polar);
                            angle = 0;
                            for a = 1:number_azimuth*2
                                for p = 1:number_polar
                                    angle = angle + 1;
                                    psi_ap(a, p) = f(:, angle)'*basis_space(:, ord_s);
                                end
                            end
                            psi_p = zeros(number_polar, 1);
                            for p = 1:number_polar
                                psi_p(p) = psi_ap(:, p)'*basis_azimuth(:, ord_a);
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

