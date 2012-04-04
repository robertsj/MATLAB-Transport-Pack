%> @file  C5G7_UO2_assembly.m
%> @brief Solve the UO2 assembly of C5G7.

clear classes;
plots = 0;

tic
% ==============================================================================
% INPUT
% ==============================================================================
input = Input();
put(input, 'number_groups',         7);
put(input, 'eigen_tolerance',       1e-4);
put(input, 'eigen_max_iters',       2);
put(input, 'outer_tolerance',       0.002);
put(input, 'outer_max_iters',       20);
put(input, 'inner_tolerance',       0.001);
put(input, 'inner_solver',          'SI');
put(input, 'livolant_free_iters',   3);
put(input, 'livolant_accel_iters',  3);
put(input, 'bc_left',               'reflect');
put(input, 'bc_right',              'reflect');
put(input, 'bc_top',                'reflect');
put(input, 'bc_bottom',             'reflect');
put(input, 'print_out',             1);


% ==============================================================================
% MATERIALS (Test two group data)
% ==============================================================================
mat = C5G7_materials(0);

% ==============================================================================
% PINS
% ==============================================================================

% Shared pin properties
pitch = 1.26; radii = 0.54; number = 9;
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

if plots
    figure(1), plot_mesh_map(assem1, 'MATERIAL')
end
%return
% ==============================================================================
% SETUP 
% ==============================================================================
state       = State(input, assem1);
quadrature  = LevelSymmetric(8);
boundary    = BoundaryMesh(input, assem1, quadrature);
q_e         = Source(assem1, 2);                  % Not initialized = not used.
q_f         = FissionSource(state, assem1, mat);  % Inititalized = used.
initialize(q_f);
toc

% ==============================================================================
% SOLVE 
% ==============================================================================
tic
solver = Eigensolver(input, state, boundary, assem1, mat, quadrature, q_e, q_f);
toc
tic
  out = solve(solver); 
toc

% ==============================================================================
% POSTPROCESS 
% ==============================================================================
%VTK.savevtk(state, 7, mesh, 'c5g7_uo2_n8_qr2.vtk')
%save('c5g7_uo2_n8_qr2.mat')
