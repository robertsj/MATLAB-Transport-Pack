%> @file  C5G7.m
%> @brief Two-dimensional, seven group transport benchmark problem.

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
mat = C5G7_materials();

% ==============================================================================
% PINS
% ==============================================================================

% Shared pin properties
pitch = 1.26; radii = 0.54; number = 8;
% Pin 1 - UO2 
matid = [1 7];  pin1 = PinCell(pitch, radii, matid); meshify(pin1, number);
% Pin 2 - 4.3% MOX
matid = [2 7];  pin2 = PinCell(pitch, radii, matid); meshify(pin2, number);
% Pin 3 - 7.0% MOX
matid = [3 7];  pin3 = PinCell(pitch, radii, matid); meshify(pin3, number);
% Pin 4 - 8.7% MOX
matid = [4 7];  pin4 = PinCell(pitch, radii, matid); meshify(pin4, number);
% Pin 5 - Guide Tube
matid = [5 7];  pin5 = PinCell(pitch, radii, matid); meshify(pin5, number);
% Pin 6 - Fission Chamber
matid = [6 7];  pin6 = PinCell(pitch, radii, matid); meshify(pin6, number);
% Pin 7 - Moderator
matid = [7  ];  pin7 = PinCell(pitch, [   ], matid); meshify(pin7, number);

% ==============================================================================
% ASSEMBLIES
% ==============================================================================

G = 5; F = 6;

% Assembly 1 - UO2
pin_map1 = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
            1 1 1 1 1 G 1 1 G 1 1 G 1 1 1 1 1 
            1 1 1 G 1 1 1 1 1 1 1 1 1 G 1 1 1 
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
            1 1 G 1 1 G 1 1 G 1 1 G 1 1 G 1 1 
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
            1 1 G 1 1 G 1 1 F 1 1 G 1 1 G 1 1 
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
            1 1 G 1 1 G 1 1 G 1 1 G 1 1 G 1 1 
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
            1 1 1 G 1 1 1 1 1 1 1 1 1 G 1 1 1 
            1 1 1 1 1 G 1 1 G 1 1 G 1 1 1 1 1 
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ];     
assem1 = Assembly({pin1, pin2, pin3, pin4, pin5, pin6, pin7}, pin_map1); 
meshify(assem1);
% Assembly 2 - MOX
pin_map2 = [2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
            2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2
            2 3 3 3 3 G 3 3 G 3 3 G 3 3 3 3 2
            2 3 3 G 3 4 4 4 4 4 4 4 3 G 3 3 2
            2 3 3 3 4 4 4 4 4 4 4 4 4 3 3 3 2
            2 3 G 4 4 G 4 4 G 4 4 G 4 4 G 3 2
            2 3 3 4 4 4 4 4 4 4 4 4 4 4 3 3 2
            2 3 3 4 4 4 4 4 4 4 4 4 4 4 3 3 2
            2 3 G 4 4 G 4 4 F 4 4 G 4 4 G 3 2
            2 3 3 4 4 4 4 4 4 4 4 4 4 4 3 3 2
            2 3 3 4 4 4 4 4 4 4 4 4 4 4 3 3 2
            2 3 G 4 4 G 4 4 G 4 4 G 4 4 G 3 2
            2 3 3 3 4 4 4 4 4 4 4 4 4 3 3 3 2
            2 3 3 G 3 4 4 4 4 4 4 4 3 G 3 3 2
            2 3 3 3 3 G 3 3 G 3 3 G 3 3 3 3 2
            2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2
            2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2];

assem2 = Assembly({pin1, pin2, pin3, pin4, pin5, pin6, pin7}, pin_map2); 
meshify(assem2);
% Assembly 3 - Moderator
pin_map3 = [7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7];
assem3 = Assembly({pin1, pin2, pin3, pin4, pin5, pin6, pin7}, pin_map3); 
meshify(assem3);
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
VTK.savevtk(state, 2, mesh, 'c5g7_n8_qr2.vtk')
save('c5g7_n8_qr2.mat')
