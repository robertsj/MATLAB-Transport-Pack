%> @file  test_fixed.m
%> @brief Tests the fixed source iteration
%
%> Finish me.
% ==============================================================================

clear

% Get the default input.
input = Input();

% Material etc.
put(input, 'number_groups',         2);

% Inner iteration parameters.
put(input, 'inner_solver',          'Livolant');
put(input, 'livolant_free_iters',   3);
put(input, 'livolant_accel_iters',  3);

% One material, one group, c = 0.9
mat         = test_materials(2);

% Simple uniformly meshed square
mesh        = test_mesh(1);

state       = State(input, mesh);
quadrature  = LevelSymmetric(8);
boundary    = Boundary(input, mesh, quadrature);

% Uniform source.
q_e         = Source(mesh, 2);
spectra     = [ 1.0
                1.0 ];      
placement   = [ 1 ];    
set_sources(q_e, spectra, placement);

% Empty fission source.
q_f = FissionSource(mesh, mat);

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
% and plot it
subplot(2,1,1)
plot_flux(mesh, f)
% Get the flux
f = flux(state, 2);
% and plot it
subplot(2,1,2)
plot_flux(mesh, f)
