%> @file  test_pin_cell.m
%> @brief Tests the eigenvalue solver on a 2D pin cell problem.
%
%> Finish me.
% ==============================================================================

clear% classes
warning('OFF', 'user:input')
% ==============================================================================
% INPUT
% ==============================================================================
input = Input();
put(input, 'number_groups',         2);
put(input, 'eigen_tolerance',       1e-5);
put(input, 'eigen_max_iters',       10);
put(input, 'inner_tolerance',       0.03);
put(input, 'inner_solver',          'SI');
put(input, 'livolant_free_iters',   3);
put(input, 'livolant_accel_iters',  3);
put(input, 'bc_left',               'reflect');
put(input, 'bc_right',              'reflect');
put(input, 'bc_top',                'reflect');
put(input, 'bc_bottom',             'reflect');

input.number_groups = 2;


% ==============================================================================
% MATERIALS (Test two group data)
% ==============================================================================
mat = test_materials(2);


% ==============================================================================
% MESH (A simple pin-cell)
% ==============================================================================

% Basic pin properties
p = 1.27;                 
r = [0.54];
% Materials: inside out
matid       = [2 1]; 
% Build the pin.
pin = PinCell(p, r, matid);
% Mesh the pin.
meshify(pin, 34);
% plot
%plot_mesh(pin);

%q = UniformMOC(5, 1, 21);
% q = CollocatedMOC(27, 1, 1);
% track(pin, q);
% clf
% Plotting tracks
%plot_tracks(pin, 0);
%verify_tracks(pin)
% 
% ==============================================================================
% SETUP 
% ==============================================================================
state       = State(input, pin);
quadrature  = LevelSymmetric(4);
boundary    = Boundary(input, pin, quadrature);
q_e         = Source(pin, 2);                  % Not initialized = not used.
q_f         = FissionSource(state, pin, mat);  % Inititalized = used.
initialize(q_f);

% ==============================================================================
% SOLVE 
% ==============================================================================
solver = Eigensolver(input,         ...
                     state,         ...
                     boundary,      ...
                     pin,          ...
                     mat,           ...
                     quadrature,	...
                     q_e,           ...
                     q_f);
 
tic
out = solve(solver); 
toc
VTK.savevtk(state, 2, pin, 'pin.vtk')
% Get the flux
f = flux(state, 1);
% and plot it
subplot(2,1,1)
plot_flux(pin, f)
axis square
% Get the flux
f = flux(state, 2);
% and plot it
subplot(2,1,2)
plot_flux(pin, f)
axis square