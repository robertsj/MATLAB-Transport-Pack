%> @file  test_response.m
%> @brief Tests a response function boundary problem.
%>
%> The problem is a two group node with multiplication (that can be turned
%> on or off).  We're looking for the outgoing fluxes.
%>
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
put(input, 'bc_left',               'reflect');
put(input, 'bc_right',              'reflect');
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

% Top

% make basis
number_space   = number_cells_x(mesh);
number_polar   = get(input, 'quad_number_polar');
number_azimuth = get(input, 'quad_number_azimuth');                
basis_space    = DiscreteLP(number_space-1);
basis_polar    = DiscreteLP(number_polar-1);
basis_azimuth  = DiscreteLP(number_azimuth*2-1);
num_ang = number_polar*number_azimuth;

octants = [2; 1];
for o = 1:2
    o_in  = octants(o, 1); % incident octant
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
    for s = 1:number_space
        angle = 1;
        for a = a1:a3:a2
            f(s, (o-1)*num_ang+angle) = ft(s, a);
            angle = angle + 1;
        end
    end
end

% Space->Azimuth->Polar.  f(space, angle, group)
coef = zeros(2*2*1,1); % first order in all.
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
            coef(i) = psi_p'*basis_polar(:, ord_p);     
            i = i + 1;
        end
    end
end



