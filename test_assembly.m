%> @file  test_assembly.m
%> @brief Tests the eigenvalue solver on a 2D assembly problem.
% 
% =========================================================================
%     Iter:     6, Error:   0.00027425, keff:   0.492221, Sweeps:   174
%                       Rate:   0.21176737
% =========================================================================
% Elapsed time is 21876.538294 seconds.
% =========================================================================
%
%
% That's about 2 minutes per sweep, wish seems reasonable enough.  However,
% if we could limit that to maybe 10 sweeps via CMFD, we're in 
% business.  Note, the GMRES

%clear

% ==============================================================================
% INPUT
% ==============================================================================
input = Input();
put(input, 'number_groups',         2);
put(input, 'number_groups',         2);
put(input, 'eigen_tolerance',       1e-4);
put(input, 'eigen_max_iters',       100);
put(input, 'inner_tolerance',       0.1);
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
% MESH (Assembly of pins)
% ==============================================================================

% Shared pin properties
pitch   = 1.26;                 
radii   = 0.54;
number  = 20;
% Pin 1
matid  = [2 1]; % IN to OUT
pin1   = PinCell(pitch, radii, matid);
meshify(pin1, number);
figure(1)
plot_mesh(pin1);
% Pin 2
matid  = [4 1]; 
pin2   = PinCell(pitch, radii, matid);
meshify(pin2, number);
figure(2)
plot_mesh(pin2);
% Pin 3
matid  = [1]; 
pin3   = PinCell(pitch, [], matid);
meshify(pin3, number);
figure(3)
plot_mesh(pin3);

% Make the assembly.
pin_map = [ 1 1 1
            1 2 1
            1 1 1];
mesh = Assembly({pin1, pin2, pin3}, pin_map);
meshify(mesh);
figure(4)
plot_mesh_map(mesh, 'MATERIAL')
return

% ==============================================================================
% SETUP 
% ==============================================================================
state       = State(input, mesh);
quadrature  = LevelSymmetric(6);
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

% ==============================================================================
% POSTPROCESS 
% ==============================================================================

figure(5)
subplot(2, 1, 1)
f1 = flux(state, 1);
plot_flux(mesh, f1)
axis equal
subplot(2, 1, 2)
f2 = flux(state, 2);
plot_flux(mesh, f2)
axis equal