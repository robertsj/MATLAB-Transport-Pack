%> @file  test_pin_cell.m
%> @brief Tests the eigenvalue solver on a 2D pin cell problem.
%
%> Finish me.
% ==============================================================================

clear
warning('OFF', 'user:input')
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
p = 1.27;                 
r = [0.54];
% Materials: inside out
matid       = [2 1 2 1]; 
% Build the pin.
pin = PinCell(p, r, matid);
% Mesh the pin.
%meshify(pin, 40);
% plot
%plot_pin(mesh);


%q = UniformMOC(5, 1, 21);
q = CollocatedMOC(27, 1, 1);
track(pin, q);
clf
% Plotting tracks
%plot_tracks(pin, 0);
verify_tracks(pin)
% 
% % ==============================================================================
% % SETUP 
% % ==============================================================================
% state       = State(input, mesh);
% quadrature  = LevelSymmetric(8);
% boundary    = Boundary(input, mesh, quadrature);
% q_e         = Source(mesh, 2);                  % Not initialized = not used.
% q_f         = FissionSource(state, mesh, mat);  % Inititalized = used.
% initialize(q_f);
% 
% % ==============================================================================
% % SOLVE 
% % ==============================================================================
% solver = Eigensolver(input,         ...
%                      state,         ...
%                      boundary,      ...
%                      mesh,          ...
%                      mat,           ...
%                      quadrature,	...
%                      q_e,           ...
%                      q_f);
%  
% tic
% out = solve(solver); 
% toc