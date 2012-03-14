%> @file  benchmark_15x15_assembly.m
%> @brief Simple 15x15 assembly.
%
%> The benchmark is a simple 15x15 consisting of a 5x5 array of 3x3 
%> subassemblies.  It's not supposed to be a very meaningful problem, but
%> it does have some gradients even when fully reflected.
%>
%> The full assembly is modeled, and a set of high order responses is 
%> generated for use with ERME.
clear;

flag = 2;

% ==============================================================================
% Quadrature (3x6 = 18 angles per quadrant)
% ==============================================================================
quadrature  = QuadrupleRange(18);

% ==============================================================================
% MATERIALS (Test two group data)
% ==============================================================================
mat = test_materials(2);

% ==============================================================================
% GEOMETRY
% ==============================================================================

% Pins
pitch = 1.26; radii = 0.54; number = 13;
% Pin 1 - Fuel 1
matid = [2 1];  pin1 = PinCell(pitch, radii, matid); meshify(pin1, number);
% Pin 2 - Fuel 2
matid = [4 1];  pin2 = PinCell(pitch, radii, matid); meshify(pin2, number);
% Pin 3 - MOD
matid = [  1];  pin3 = PinCell(pitch,    [], matid); meshify(pin3, number);

% Assemblies
pin_map1 = [1 1 1 
            1 1 1 
            1 1 1];     
assem1 = Assembly({pin1, pin2, pin3}, pin_map1); % kinf = 1.334776  
meshify(assem1);
pin_map2 = [1 1 1 
            1 2 1 
            1 1 1];     
assem2 = Assembly({pin1, pin2, pin3}, pin_map2); % kinf = 0.906918
meshify(assem2);
pin_map3 = [1 1 1 
            1 3 1 
            1 1 1];     
assem3 = Assembly({pin1, pin2, pin3}, pin_map3); % kinf = 1.319136
meshify(assem3);

% ==============================================================================
% Benchmark Solve
% ==============================================================================
if flag == 0

input = Input();
put(input, 'number_groups',         2);
put(input, 'eigen_tolerance',       1e-8);
put(input, 'eigen_max_iters',       30);
put(input, 'inner_tolerance',       1e-9);
put(input, 'inner_max_iters',       100);
put(input, 'outer_tolerance',       1e-9);
put(input, 'outer_max_iters',       20);
put(input, 'inner_solver',          'SI');
put(input, 'bc_left',               'reflect');
put(input, 'bc_right',              'reflect');
put(input, 'bc_top',                'reflect');
put(input, 'bc_bottom',             'reflect');

% Core
assembly_map = [ 3 2 1 2 1
                 2 1 2 1 2
                 1 2 1 2 1
                 2 1 2 1 2
                 1 2 1 2 1 ];
core = Core({assem1, assem2, assem3}, assembly_map); 
meshify(core);
if 1==0, figure(1), plot_mesh_map(core, 'MATERIAL'), shading flat, end

state       = State(input, core);
boundary    = BoundaryMesh(input, core, quadrature);
q_e         = Source(core, 2);
q_f = FissionSource(state, core, mat);
initialize(q_f);
solver= Eigensolver(input,      ...
                    state,    	...
                    boundary,   ...
                    core,     	...
                    mat,        ...
                    quadrature, ...
                    q_e,        ...
                    q_f);
solve(solver);          

f1 = flux(state, 1); figure(1), plot_flux(core, f1), shading flat;
f2 = flux(state, 2); figure(2), plot_flux(core, f2), shading flat;
VTK.savevtk(state, 2, core, 'benchmark_15x15.vtk');
save('benchmark_15x15.mat')
% keff = 1.103240559262379

% ==============================================================================
% Generate DB
% ==============================================================================

elseif flag == 1 

input = Input();
put(input, 'number_groups',         2);
put(input, 'inner_tolerance',       1e-12);
put(input, 'inner_max_iters',       100);
put(input, 'outer_tolerance',       1e-12);
put(input, 'outer_max_iters',       10);
put(input, 'rf_max_order_space',    4);
put(input, 'rf_max_order_azimuth',  4); % max is 11
put(input, 'rf_max_order_polar',    2); % max is 2
put(input, 'quad_order',           18); 
put(input, 'rf_k_vector',           [0.7 0.8 0.9 1.0 1.1 1.2]');
put(input, 'rf_number_nodes',       3);
put(input, 'rf_db_name',            'benchmark_15x15_rf_db.h5');

% Make and run a driver
mesh_array = {assem1, assem2, assem2};
driver = ResponseDriver(input, mat, mesh_array);
run(driver);
[R, F, A, L] = get_responses(driver);

% Make and write a database
rf_db = ResponseDB(input);
initialize_write(rf_db);
write_response(rf_db, 1, 'fuel_3x3',            ...
    R(:,:,:,1), F(:,:,1), A(:,:,1), L(:,:,:,1));
write_response(rf_db, 2, 'fuel_gd_3x3',         ...
    R(:,:,:,2), F(:,:,2), A(:,:,2), L(:,:,:,2));
write_response(rf_db, 3, 'fuel_water_hole_3x3', ...
    R(:,:,:,3), F(:,:,3), A(:,:,3), L(:,:,:,3));


% ==============================================================================
% RUN THE ERME PROBLEM
% ==============================================================================
else

disp(' doing benchmark erme ')
    
% ==============================================================================
% Server (here, a the DB)
% ==============================================================================
db_input = Input();
put(db_input, 'rf_db_name', 'benchmark_15x15_rf_db.h5');
rf_db = ResponseDB(db_input);
read_response(rf_db);

% ==============================================================================
% ERME setup
% ==============================================================================
input = Input();
put(input, 'bc_left',           'reflect');
put(input, 'bc_right',          'reflect');
put(input, 'bc_bottom',         'reflect');
put(input, 'bc_top',            'reflect');
put(input, 'inner_tolerance',   1e-12);
put(input, 'inner_max_iters',   100);
put(input, 'outer_tolerance',   1e-7);
put(input, 'outer_max_iters',   10);
put(input, 'rf_order_space',    4); % Eventually, these will be the
put(input, 'rf_order_azimuth',  4); % orders extracted from a DB
put(input, 'rf_order_polar',    2); 
put(input, 'number_groups',     2);
put(input, 'dimension',         2);

assembly_map = [ 3 2 1 2 1
                 2 1 2 1 2
                 1 2 1 2 1
                 2 1 2 1 2
                 1 2 1 2 1 ];
%assembly_map = [ 2 ];
problem = ERME(input, rf_db, assembly_map);
solver = ERME_Picard(input, problem);
solve(solver);

end


