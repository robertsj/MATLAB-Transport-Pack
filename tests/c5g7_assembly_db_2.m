%> @file  c5g7_pin_db_erme.m
%> @brief Generate response database for C5G7 assemblies cells.
%

clear
tic
% Include the MTP source directory
path(path, '../')
Mesh();
% ==============================================================================
% Create a response DB
% ==============================================================================

inputrf = Input();
put(inputrf, 'number_groups',         7);
put(inputrf, 'dimension',             2);
put(inputrf, 'inner_tolerance',       1e-12);
put(inputrf, 'inner_max_iters',       100);
put(inputrf, 'outer_tolerance',       1e-9);
put(inputrf, 'outer_max_iters',       10);
put(inputrf, 'inner_print_out',       1);
put(inputrf, 'rf_print_out',          1);
put(inputrf, 'quad_order',           18);
put(inputrf,  'pc_switch',              1);
put(inputrf, 'outer_precondition',    1);
put(inputrf, 'rf_max_order_space',    4);
put(inputrf, 'rf_max_order_azimuth',  4);
put(inputrf, 'rf_max_order_polar',    2);
% Different keff vector for the moderator (i.e. 1 value)
keff_vectors    = cell(1, 1);
keff_vectors{1} = [0.9; 1.0; 1.1; 1.2; 1.3];
put(inputrf, 'rf_keff_vector',        keff_vectors);
put(inputrf, 'rf_number_nodes',       1);
put(inputrf, 'rf_db_name',            'c5g7_assemblies2.h5');
descriptions = {'MOX'};
put(inputrf, 'rf_node_descriptions',  descriptions);
% ==============================================================================
% MATERIALS
% ==============================================================================
mat = C5G7_materials(0);

% ==============================================================================
% PINS
% ==============================================================================

% Shared pin properties
pitch = 1.26; radii = 0.54; number = 13; %52;
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

% ==============================================================================
% RESPONSE FUNCTION DRIVER
% ==============================================================================
mesh_array = {assem2};


driver = ResponseDriver(inputrf, mat, mesh_array);
run(driver);





