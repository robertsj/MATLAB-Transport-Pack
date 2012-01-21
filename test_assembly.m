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
matid  = [1 2]; 
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
% Pin 3
matid  = [1]; 
pin3   = PinCell(pitch, [], matid);
meshify(pin3, number);
figure(3)
plot_pin(pin3);

% Make the assembly.
pin_map = [ 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3     
            3 3 3 1 2 1 2 1 2 1 2 1 2 1 3 3 3 
            3 3 1 2 1 2 1 2 1 2 1 2 1 2 1 3 3   
            3 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 3 
            3 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 3   
            3 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 3 
            3 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 3   
            3 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 3 
            3 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 3   
            3 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 3 
            3 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 3   
            3 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 3 
            3 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 3   
            3 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 3 
            3 3 1 2 1 2 1 2 1 2 1 2 1 2 1 3 3   
            3 3 3 1 2 1 2 1 2 1 2 1 2 1 3 3 3 
            3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 ];
mesh = Assembly({pin1, pin2, pin3}, pin_map);
meshify(mesh);
figure(4)
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