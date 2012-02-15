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
put(input, 'inner_solver',          'SI');
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
put(input, 'quad_number_polar',     2);
put(input, 'quad_number_azimuth',   4);
% One material, one group, c = 0.9
mat         = test_materials(1);

% Simple uniformly meshed square
mesh        = test_mesh(1);

state       = State(input, mesh);
quadrature  = QuadrupleRange(8);
boundary    = BoundaryMesh(input, mesh, quadrature);

% Uniform source.
q_e         = Source(mesh, 2);
spectra     = [ 1.0 0
                1.0 0];      
placement   = [ 2 2; 2 2];    
set_sources(q_e, spectra, placement);

% Empty fission source.
q_f = FissionSource(state, mesh, mat);
%initialize(q_f);

% Make the inner iteration.
solver= Fixed(input,        ...
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

% Get the flux
f = flux(state, 1);
subplot(2,1,1)
plot_flux(mesh, f)
axis square
% Get the flux
f = flux(state, 1);
% and plot it
subplot(2,1,2)
plot_flux(mesh, f)
axis square
