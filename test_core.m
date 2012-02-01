%> @file  test_core.m
%> @brief Tests the eigenvalue solver on a 2D core.

clear all;
plots = 0;

% ==============================================================================
% INPUT
% ==============================================================================
input = Input();
put(input, 'number_groups',         2);
put(input, 'number_groups',         2);
put(input, 'eigen_tolerance',       1e-4);
put(input, 'eigen_max_iters',       100);
put(input, 'inner_tolerance',       0.0001);
put(input, 'inner_solver',          'SI');
put(input, 'livolant_free_iters',   3);
put(input, 'livolant_accel_iters',  3);
put(input, 'bc_left',               'reflect');
put(input, 'bc_right',              'vacuum');
put(input, 'bc_top',                'reflect');
put(input, 'bc_bottom',             'vacuum');
input.number_groups = 2;

% ==============================================================================
% MATERIALS (Test two group data)
% ==============================================================================
mat = test_materials(2);

% ==============================================================================
% PINS
% ==============================================================================

% Shared pin properties
pitch = 1.26; radii = 0.54; number = 8;
% Pin 1
matid = [2 1];  pin1 = PinCell(pitch, radii, matid); meshify(pin1, number);
% Pin 2
matid = [4 1];  pin2 = PinCell(pitch, radii, matid); meshify(pin2, number);
% Pin 3
matid = [1];    pin3 = PinCell(pitch, [   ], matid); meshify(pin3, number);
if plots
    figure(1); plot_mesh(pin1);
    figure(2), plot_mesh(pin2);
    figure(3), plot_mesh(pin3);
end

% ==============================================================================
% ASSEMBLIES
% ==============================================================================

% Assembly 1
pin_map1 = [1 1 1; 1 2 1; 1 1 1];
assem1 = Assembly({pin1, pin2, pin3}, pin_map1); meshify(assem1);
% Assembly 2
pin_map2 = [1 1 1; 1 1 1; 1 1 1];
assem2 = Assembly({pin1, pin2, pin3}, pin_map2); meshify(assem2);
% Assembly 3
pin_map3 = [3 3 3; 3 3 3; 3 3 3];
assem3 = Assembly({pin1, pin2, pin3}, pin_map3); meshify(assem3);
if plots
    figure(4), plot_mesh_map(assem1, 'MATERIAL')
    figure(5), plot_mesh_map(assem2, 'MATERIAL')
    figure(6), plot_mesh_map(assem3, 'MATERIAL')
end

% ==============================================================================
% CORE
% ==============================================================================

assembly_map = [ 1 2 3; 2 1 3; 3 3 3];
core = Assembly({assem1, assem2, assem3}, assembly_map); meshify(core);
if plots, figure(7), plot_mesh_map(core, 'MATERIAL'), end


% ==============================================================================
% SETUP 
% ==============================================================================
state       = State(input, core);
quadrature  = QuadrupleRange(2);
boundary    = BoundaryMesh(input, core, quadrature);
q_e         = Source(core, 2);                  % Not initialized = not used.
q_f         = FissionSource(state, core, mat);  % Inititalized = used.
initialize(q_f);

% ==============================================================================
% SOLVE 
% ==============================================================================
solver = Eigensolver(input, state, boundary, core, mat, quadrature, q_e, q_f);
tic
  out = solve(solver); 
toc

% ==============================================================================
% POSTPROCESS 
% ==============================================================================
VTK.savevtk(state, 2, mesh, 'test_core.vtk')
save('test_core.mat')
