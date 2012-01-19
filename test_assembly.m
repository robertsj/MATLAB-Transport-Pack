%> @file  test_assembly.m
%> @brief Tests the eigenvalue solver on a 2D assembly problem.
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
% MESH (Assembly of pins)
% ==============================================================================

% Shared pin properties
pitch   = 1.26;                 
radii   = 0.54;
number  = 14;
% Pin 1
matid  = [4 2]; 
pin1   = PinCell(pitch, radii, matid);
meshify(pin1, number);
figure(1)
plot_pin(pin1);
% Pin 2
matid  = [1 3]; 
pin2   = PinCell(pitch, radii, matid);
meshify(pin2, number);
figure(2)
plot_pin(pin2);

% Make the assembly.
pin_map = [ 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1   
            2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 
            1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1   
            2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 
            1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1   
            2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 
            1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1   
            2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 
            1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1   
            2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 
            1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1   
            2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 
            1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1   
            2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 
            1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1   
            2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 
            1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 ];
mesh = Assembly({pin1, pin2}, pin_map);
meshify(mesh);
figure(3)
plot_mesh_map(mesh, 'MATERIAL')


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

% ==============================================================================
% POSTPROCESS 
% ==============================================================================

figure(4)
subplot(2, 1, 1)
f1 = flux(state, 1);
plot_flux(mesh, f1)
axis equal
subplot(2, 1, 2)
f2 = flux(state, 2);
plot_flux(mesh, f2)
axis equal