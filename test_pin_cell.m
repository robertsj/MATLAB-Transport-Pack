%> @file  test_pin_cell.m
%> @brief Tests the eigenvalue solver on a 2D pin cell problem.
%
%> Finish me.
% ==============================================================================

clear

% ==============================================================================
% INPUT
% ==============================================================================
input = Input();
put(input, 'number_groups',         2);
put(input, 'inner_solver',          'GMRES');
put(input, 'livolant_free_iters',   3);
put(input, 'livolant_accel_iters',  3);
input.number_groups = 2;


% ==============================================================================
% MATERIALS (Test two group data)
% ==============================================================================
mat = test_materials(2);


% ==============================================================================
% MESH (A simple pin-cell)
% ==============================================================================

% Basic pin properties
pitch       = 20.0;                 
radii       = [2 5 8]
% Materials: inside out
matid       = [2 1 2 1]; 
% Build the pin.
mesh = PinCell(pitch, radii, matid);
% Mesh the pin.
meshify(mesh, 40);
% plot
plot_pin(mesh);

% ==============================================================================
% SETUP 
% ==============================================================================
state       = State(input, mesh);
quadrature  = LevelSymmetric(8);
boundary    = Boundary(input, mesh, quadrature);
q_e         = Source(mesh, 2);                  % Not initialized = not used.
q_f         = FissionSource(state, mesh, mat);  % Inititalized = used.
initialize(q_f);

% ==============================================================================
% SOLVE 
% ==============================================================================
solver = Eigensolver(input,         ...
                     state,         ...
                     boundary,      ...
                     mesh,          ...
                     mat,           ...
                     quadrature,	...
                     q_e,           ...
                     q_f);
 
tic
out = solve(solver); 
toc